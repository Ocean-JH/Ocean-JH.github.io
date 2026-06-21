#!/usr/bin/env python3
"""Manage repeatable website content without external dependencies."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTENT_PATH = ROOT / "assets/data/content.json"
POST_CONTENT_PATH = ROOT / "assets/data/post-content.json"
NEWS_CONTENT_PATH = ROOT / "assets/data/news-content.json"
POSTS_DIR = ROOT / "content/posts"
NEWS_DIR = ROOT / "content/news"


def read_json(path: Path, default):
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def json_text(value) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2) + "\n"


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json_text(value), encoding="utf-8")


def load_content() -> dict:
    if CONTENT_PATH.exists():
        data = read_json(CONTENT_PATH, {})
    else:
        data = {
            "schemaVersion": 1,
            "news": [],
            "research": [],
            "publications": [],
            "posts": [],
        }

    data.setdefault("schemaVersion", 1)
    for key in ("news", "research", "publications", "posts"):
        data.setdefault(key, [])
        if not isinstance(data[key], list):
            raise ValueError(f"content.{key} must be a list")
    return data


def parse_scalar(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        return str(json.loads(value))
    if len(value) >= 2 and value[0] == value[-1] == "'":
        return value[1:-1].replace("''", "'")
    return value


def parse_front_matter(raw: str, path: Path) -> tuple[dict, str]:
    if not raw.startswith("---"):
        return {}, raw

    opening_end = raw.find("\n")
    closing_start = raw.find("\n---", opening_end + 1)
    if opening_end < 0 or closing_start < 0:
        raise ValueError(f"Malformed front matter in {path.relative_to(ROOT)}")

    header = raw[opening_end + 1 : closing_start]
    body_start = closing_start + 4
    if raw[body_start : body_start + 2] == "\r\n":
        body_start += 2
    elif raw[body_start : body_start + 1] == "\n":
        body_start += 1

    metadata: dict[str, object] = {}
    active_list = ""
    for line in header.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if active_list and stripped.startswith("- "):
            metadata[active_list].append(parse_scalar(stripped[2:]))
            continue

        match = re.match(r"^([A-Za-z][A-Za-z0-9_-]*):\s*(.*)$", stripped)
        if not match:
            continue
        key, value = match.groups()
        if not value:
            metadata[key] = [] if key == "tags" else ""
            active_list = key if key == "tags" else ""
        elif key == "tags" and value.startswith("["):
            metadata[key] = [parse_scalar(item) for item in value[1:-1].split(",") if item.strip()]
            active_list = ""
        else:
            metadata[key] = parse_scalar(value)
            active_list = ""

    return metadata, raw[body_start:]


def normalize_date(value: str) -> str:
    match = re.fullmatch(r"(\d{4})-(\d{1,2})-(\d{1,2})", str(value).strip())
    if not match:
        raise ValueError(f"Invalid date {value!r}; use YYYY-MM-DD")
    year, month, day = map(int, match.groups())
    return f"{year:04d}-{month:02d}-{day:02d}"


def slug_from_path(path: Path) -> str:
    return re.sub(r"^\d{4}-\d{2}-\d{2}-", "", path.stem)


def build_posts(data: dict) -> tuple[list[dict], dict[str, str]]:
    existing = {str(item.get("slug", "")): item for item in data.get("posts", [])}
    posts: list[dict] = []
    bodies: dict[str, str] = {}
    seen: set[str] = set()

    for path in sorted(POSTS_DIR.glob("*.md")):
        raw = path.read_text(encoding="utf-8")
        front_matter, _ = parse_front_matter(raw, path)
        filename_slug = slug_from_path(path)
        slug = str(front_matter.get("slug") or filename_slug).strip()
        fallback = existing.get(slug, {})

        if not slug or slug in seen:
            raise ValueError(f"Missing or duplicate post slug: {slug!r}")
        seen.add(slug)

        filename_date = path.name[:10]
        date = normalize_date(str(front_matter.get("date") or fallback.get("date") or filename_date))
        tags = front_matter.get("tags") or fallback.get("tags") or []
        if not isinstance(tags, list):
            raise ValueError(f"tags must be a list in {path.relative_to(ROOT)}")

        post = {
            "title": str(front_matter.get("title") or fallback.get("title") or "").strip(),
            "slug": slug,
            "date": date,
            "topic": str(front_matter.get("topic") or fallback.get("topic") or "").strip(),
            "tags": [str(tag).strip() for tag in tags if str(tag).strip()],
            "description": str(front_matter.get("description") or fallback.get("description") or "").strip(),
            "source": f"/content/posts/{path.name}",
            "url": f"/posts/post.html?slug={slug}",
        }

        missing = [key for key in ("title", "date", "topic", "description") if not post[key]]
        if missing:
            fields = ", ".join(missing)
            raise ValueError(f"Missing {fields} in {path.relative_to(ROOT)}")

        posts.append(post)
        bodies[slug] = raw

    posts.sort(key=lambda item: (item["date"], item["title"].lower()), reverse=True)
    return posts, bodies


def build_news(data: dict) -> tuple[list[dict], dict[str, str]]:
    existing_items = data.get("news", [])
    existing = {str(item.get("slug", "")): item for item in existing_items if item.get("slug")}
    news = [item for item in existing_items if not item.get("source")]
    bodies: dict[str, str] = {}
    seen = {str(item.get("slug", "")) for item in news if item.get("slug")}

    for path in sorted(NEWS_DIR.glob("*.md")):
        raw = path.read_text(encoding="utf-8")
        front_matter, _ = parse_front_matter(raw, path)
        slug = str(front_matter.get("slug") or slug_from_path(path)).strip()
        fallback = existing.get(slug, {})

        if not slug or slug in seen:
            raise ValueError(f"Missing or duplicate news slug: {slug!r}")
        seen.add(slug)

        filename_date = path.name[:10]
        item = {
            "title": str(front_matter.get("title") or fallback.get("title") or "").strip(),
            "slug": slug,
            "date": normalize_date(str(front_matter.get("date") or fallback.get("date") or filename_date)),
            "description": str(front_matter.get("description") or fallback.get("description") or "").strip(),
            "image": str(front_matter.get("image") or fallback.get("image") or "").strip(),
            "imageAlt": str(front_matter.get("imageAlt") or fallback.get("imageAlt") or "").strip(),
            "source": f"/content/news/{path.name}",
            "url": f"/news/news.html?slug={slug}",
        }
        news.append(item)
        bodies[slug] = raw

    news.sort(key=lambda item: (item.get("date", ""), item.get("title", "").lower()), reverse=True)
    return news, bodies


def validate_collections(data: dict) -> None:
    requirements = {
        "news": ("date", "title", "image", "imageAlt", "url"),
        "research": ("title", "description"),
        "publications": ("title", "authors", "venue", "year"),
    }
    for collection, required_fields in requirements.items():
        for index, item in enumerate(data.get(collection, [])):
            if not isinstance(item, dict):
                raise ValueError(f"content.{collection}[{index}] must be an object")
            missing = [field for field in required_fields if not str(item.get(field, "")).strip()]
            if missing:
                fields = ", ".join(missing)
                raise ValueError(f"content.{collection}[{index}] is missing {fields}")
            if collection == "news":
                normalize_date(str(item["date"]))
                image = str(item["image"]).strip()
                if image.startswith("/") and not image.startswith("//"):
                    image_path = ROOT / image.lstrip("/")
                    if not image_path.is_file():
                        raise ValueError(f"News image does not exist: {image}")


def generated_content(data: dict) -> tuple[dict, dict[str, str], dict[str, str]]:
    news, news_bodies = build_news(data)
    posts, bodies = build_posts(data)
    generated = {
        "schemaVersion": 1,
        "news": news,
        "research": data.get("research", []),
        "publications": data.get("publications", []),
        "posts": posts,
    }
    validate_collections(generated)
    return generated, bodies, news_bodies


def sync(check: bool = False) -> None:
    data, bodies, news_bodies = generated_content(load_content())
    if check:
        stale = []
        if read_json(CONTENT_PATH, {}) != data:
            stale.append(CONTENT_PATH.relative_to(ROOT))
        if read_json(POST_CONTENT_PATH, {}) != bodies:
            stale.append(POST_CONTENT_PATH.relative_to(ROOT))
        if read_json(NEWS_CONTENT_PATH, {}) != news_bodies:
            stale.append(NEWS_CONTENT_PATH.relative_to(ROOT))
        if stale:
            raise ValueError("Run `python3 scripts/content.py sync`; stale: " + ", ".join(map(str, stale)))
        print("Content artifacts are current.")
        return

    write_json(CONTENT_PATH, data)
    write_json(POST_CONTENT_PATH, bodies)
    write_json(NEWS_CONTENT_PATH, news_bodies)
    print(f"Synced {len(data['posts'])} posts and {len(data['news'])} news entries.")


def yaml_scalar(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def new_post(args) -> None:
    date = normalize_date(args.date)
    slug = re.sub(r"[^a-z0-9]+", "-", args.slug.lower()).strip("-")
    if not slug:
        raise ValueError("Post slug cannot be empty")
    path = POSTS_DIR / f"{date}-{slug}.md"
    if path.exists():
        raise ValueError(f"Post already exists: {path.relative_to(ROOT)}")

    tags = "\n".join(f"  - {yaml_scalar(tag)}" for tag in args.tag)
    front_matter = [
        "---",
        f"title: {yaml_scalar(args.title)}",
        f"date: {date}",
        f"topic: {yaml_scalar(args.topic)}",
        f"description: {yaml_scalar(args.description)}",
        "tags:",
        tags,
        "---",
        "",
        f"# {args.title}",
        "",
    ]
    POSTS_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(front_matter), encoding="utf-8")
    sync()
    print(f"Created {path.relative_to(ROOT)}")


def new_news(args) -> None:
    date = normalize_date(args.date)
    slug = re.sub(r"[^a-z0-9]+", "-", args.slug.lower()).strip("-")
    if not slug:
        raise ValueError("News slug cannot be empty")
    path = NEWS_DIR / f"{date}-{slug}.md"
    if path.exists():
        raise ValueError(f"News entry already exists: {path.relative_to(ROOT)}")

    front_matter = [
        "---",
        f"title: {yaml_scalar(args.title)}",
        f"date: {date}",
        f"description: {yaml_scalar(args.description)}",
        f"image: {yaml_scalar(args.image)}",
        f"imageAlt: {yaml_scalar(args.image_alt)}",
        "---",
        "",
        f"# {args.title}",
        "",
    ]
    NEWS_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(front_matter), encoding="utf-8")
    sync()
    print(f"Created {path.relative_to(ROOT)}")


def add_news(args) -> None:
    data = load_content()
    data["news"].append({
        "date": normalize_date(args.date),
        "title": args.title,
        "description": args.description,
        "image": args.image,
        "imageAlt": args.image_alt,
        "url": args.url,
    })
    data["news"].sort(key=lambda item: item.get("date", ""), reverse=True)
    write_json(CONTENT_PATH, data)
    print(f"Added news: {args.title}")


def add_research(args) -> None:
    data = load_content()
    data["research"].append({
        "title": args.title,
        "description": args.description,
        "url": args.url,
    })
    write_json(CONTENT_PATH, data)
    print(f"Added research item: {args.title}")


def add_publication(args) -> None:
    data = load_content()
    data["publications"].append({
        "title": args.title,
        "authors": args.authors,
        "venue": args.venue,
        "year": args.year,
        "doi": args.doi,
        "url": args.url,
        "pdf": args.pdf,
        "selected": args.selected,
    })
    write_json(CONTENT_PATH, data)
    print(f"Added publication: {args.title}")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)
    commands.add_parser("sync", help="Rebuild post metadata and embedded Markdown")
    commands.add_parser("check", help="Verify generated content artifacts")

    post = commands.add_parser("new-post", help="Create and index a Markdown post")
    post.add_argument("--date", required=True, help="YYYY-MM-DD")
    post.add_argument("--slug", required=True)
    post.add_argument("--title", required=True)
    post.add_argument("--topic", required=True)
    post.add_argument("--description", required=True)
    post.add_argument("--tag", action="append", default=[])

    local_news = commands.add_parser("new-news", help="Create a Markdown news entry with a detail page")
    local_news.add_argument("--date", required=True, help="YYYY-MM-DD")
    local_news.add_argument("--slug", required=True)
    local_news.add_argument("--title", required=True)
    local_news.add_argument("--description", default="")
    local_news.add_argument("--image", required=True)
    local_news.add_argument("--image-alt", required=True)

    news = commands.add_parser("add-news", help="Add news linked to an existing detail page")
    news.add_argument("--date", required=True, help="YYYY-MM-DD")
    news.add_argument("--title", required=True)
    news.add_argument("--description", default="")
    news.add_argument("--image", required=True)
    news.add_argument("--image-alt", required=True)
    news.add_argument("--url", required=True)

    research = commands.add_parser("add-research", help="Add a research item")
    research.add_argument("--title", required=True)
    research.add_argument("--description", required=True)
    research.add_argument("--url", default="")

    publication = commands.add_parser("add-publication", help="Add a publication")
    publication.add_argument("--title", required=True)
    publication.add_argument("--authors", required=True)
    publication.add_argument("--venue", required=True)
    publication.add_argument("--year", required=True)
    publication.add_argument("--doi", default="")
    publication.add_argument("--url", default="")
    publication.add_argument("--pdf", default="")
    publication.add_argument("--selected", action="store_true")
    return root


def main() -> int:
    args = parser().parse_args()
    actions = {
        "sync": lambda _: sync(),
        "check": lambda _: sync(check=True),
        "new-post": new_post,
        "new-news": new_news,
        "add-news": add_news,
        "add-research": add_research,
        "add-publication": add_publication,
    }
    try:
        actions[args.command](args)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
