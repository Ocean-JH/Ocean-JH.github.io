# Personal Academic Website

This is a no-build static site for GitHub Pages. The browser reads committed JSON and Markdown-derived artifacts directly; deployment does not run Node, Ruby, or Python.

## Content Organization

Use `assets/data/site.json` for stable site configuration:

- site metadata and navigation labels
- profile and social links
- homepage copy and section headings
- resource links
- archive and post-reader labels

Use `assets/data/content.json` for repeatable, frequently updated collections:

- `news`
- `research`
- `publications`
- generated `posts` metadata

Every collection item uses `title`, `description`, and optional `url` where applicable. Section-specific fields include `date`, publication details, post tags, and post source paths. Empty collections hide their corresponding sections.

## Content Commands

The helper uses only the Python standard library.

```bash
# Verify generated post and news artifacts
python3 scripts/content.py check

# Rebuild metadata and embedded Markdown after editing content
python3 scripts/content.py sync
```

Create a post and update both generated artifacts in one command:

```bash
python3 scripts/content.py new-post \
  --date 2026-06-20 \
  --slug example-update \
  --title "Example update" \
  --topic "Research" \
  --description "A short archive summary." \
  --tag "AI for Materials"
```

Publish news with a homepage image and local detail page:

```bash
python3 scripts/content.py new-news \
  --date 2026-06-21 \
  --slug conference-presentation \
  --title "Conference presentation" \
  --description "Presented our recent work." \
  --image "/assets/images/news/2026-06-21-conference.webp" \
  --image-alt "Wang Jianghai presenting at the conference"
```

Place the image in `assets/images/news/` first. The command creates `content/news/YYYY-MM-DD-slug.md`; edit that Markdown file with the detailed content and run `python3 scripts/content.py sync` again.

Use `add-news` when the detailed content already exists elsewhere:

```bash
python3 scripts/content.py add-news \
  --date 2026-06-21 \
  --title "External update" \
  --description "Optional detail." \
  --image "/assets/images/news/2026-06-21-update.webp" \
  --image-alt "Description of the news image" \
  --url "https://example.com/details"
```

Add research and publication entries without hand-editing JSON:

```bash
python3 scripts/content.py add-research \
  --title "Research direction" \
  --description "Short description." \
  --url ""

python3 scripts/content.py add-publication \
  --title "Paper title" \
  --authors "Author list" \
  --venue "Journal" \
  --year 2026 \
  --doi "" \
  --url "" \
  --pdf "" \
  --selected
```

## Post Files

- `content/posts/`: editable Markdown source files
- `assets/data/content.json`: collection data and generated post metadata
- `assets/data/post-content.json`: generated browser-ready Markdown map

Commit both generated JSON files after adding or editing a post. Run `python3 scripts/content.py check` before publishing.

## News Files

- `content/news/`: editable Markdown detail pages
- `assets/images/news/`: homepage and detail-page images
- `assets/data/content.json`: generated news metadata and URLs
- `assets/data/news-content.json`: generated browser-ready Markdown map

Commit the Markdown, image, `content.json`, and `news-content.json` after publishing news.
