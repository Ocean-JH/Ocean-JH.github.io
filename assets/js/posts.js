const POSTS_INDEX_URL = "/assets/data/posts.json";
const POSTS_CONTENT_URL = "/assets/data/post-content.json";
const MARKED_URL = "https://cdn.jsdelivr.net/npm/marked/marked.min.js";
const MATHJAX_URL = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js";

const cleanText = (value) => (value || "").toString().trim();

function escapeHtml(value) {
  return cleanText(value)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function normalizeSlugFromFilename(filename) {
  return cleanText(filename)
    .replace(/^.*\//, "")
    .replace(/^\d{4}-\d{1,2}-\d{1,2}-/, "")
    .replace(/\.md$/, "");
}

function formatDate(value) {
  if (!value) {
    return "";
  }
  const [year, month, day] = cleanText(value).split("-").map(Number);
  if (!year || !month || !day) {
    return cleanText(value);
  }
  return new Intl.DateTimeFormat("en", {
    year: "numeric",
    month: "short",
    day: "numeric"
  }).format(new Date(Date.UTC(year, month - 1, day)));
}

function getPostSlug() {
  return new URLSearchParams(window.location.search).get("slug") || "";
}

async function fetchJson(url) {
  const response = await fetch(url, { cache: "no-store" });
  if (!response.ok) {
    throw new Error(`Unable to load ${url}`);
  }
  return response.json();
}

async function fetchText(url) {
  const response = await fetch(url, { cache: "no-store" });
  if (!response.ok) {
    throw new Error(`Unable to load ${url}`);
  }
  return response.text();
}

function loadScript(url, timeout = 800) {
  if (document.querySelector(`script[src="${url}"]`)) {
    return Promise.resolve();
  }

  return new Promise((resolve, reject) => {
    const script = document.createElement("script");
    const timer = window.setTimeout(() => reject(new Error(`Timed out loading ${url}`)), timeout);
    script.src = url;
    script.async = true;
    script.onload = () => {
      window.clearTimeout(timer);
      resolve();
    };
    script.onerror = () => {
      window.clearTimeout(timer);
      reject(new Error(`Unable to load ${url}`));
    };
    document.head.append(script);
  });
}

async function ensureMarked() {
  if (window.marked && typeof window.marked.parse === "function") {
    return true;
  }
  try {
    await loadScript(MARKED_URL);
  } catch (error) {
    console.warn(error.message);
  }
  return Boolean(window.marked && typeof window.marked.parse === "function");
}

async function ensureMathJax() {
  window.MathJax = window.MathJax || {
    tex: {
      inlineMath: [["$", "$"], ["\\(", "\\)"]],
      displayMath: [["$$", "$$"], ["\\[", "\\]"]]
    },
    svg: { fontCache: "global" }
  };

  if (window.MathJax && typeof window.MathJax.typesetPromise === "function") {
    return true;
  }
  try {
    await loadScript(MATHJAX_URL);
  } catch (error) {
    console.warn(error.message);
  }
  return Boolean(window.MathJax && typeof window.MathJax.typesetPromise === "function");
}

function rewriteLegacyMarkdown(markdown) {
  return markdown
    .replace(/\]\(\/_posts\/\d{4}-\d{1,2}-\d{1,2}-([^)]+?)\.md\)/g, "](/posts/post.html?slug=$1)")
    .replace(/\]\(_posts\/\d{4}-\d{1,2}-\d{1,2}-([^)]+?)\.md\)/g, "](/posts/post.html?slug=$1)")
    .replace(/\]\(\.\.\/images\//g, "](/images/");
}

function stripFrontMatter(markdown) {
  return markdown.replace(/^---[\s\S]*?---\s*/, "");
}

function normalizeHeadingText(value) {
  return cleanText(value)
    .replace(/!\[[^\]]*\]\([^)]*\)/g, " ")
    .replace(/\[([^\]]+)\]\([^)]*\)/g, "$1")
    .replace(/<[^>]+>/g, " ")
    .replace(/[`*_~]/g, "")
    .replace(/\s+/g, " ")
    .trim()
    .toLowerCase();
}

function stripDuplicateTitleHeading(markdown, title) {
  const expectedTitle = normalizeHeadingText(title);
  if (!expectedTitle) {
    return markdown;
  }

  const lines = markdown.split(/\r?\n/);
  const firstContentIndex = lines.findIndex((line) => line.trim());
  if (firstContentIndex === -1) {
    return markdown;
  }

  const firstLine = lines[firstContentIndex].trim();
  const heading = /^#\s+(.+?)\s*#*$/.exec(firstLine);
  if (!heading || normalizeHeadingText(heading[1]) !== expectedTitle) {
    return markdown;
  }

  lines.splice(firstContentIndex, 1);
  return lines.join("\n").replace(/^\s+/, "");
}

function inlineMarkdown(value) {
  let output = escapeHtml(value);
  output = output.replace(/!\[([^\]]*)\]\(([^)\s]+)(?:\s+"[^"]*")?\)/g, '<img src="$2" alt="$1" loading="lazy">');
  output = output.replace(/\[([^\]]+)\]\(([^)\s]+)(?:\s+"[^"]*")?\)/g, '<a href="$2" rel="noopener">$1</a>');
  output = output.replace(/`([^`]+)`/g, "<code>$1</code>");
  output = output.replace(/\*\*\*([^*]+)\*\*\*/g, "<strong><em>$1</em></strong>");
  output = output.replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>");
  output = output.replace(/\*([^*]+)\*/g, "<em>$1</em>");
  return output;
}

function renderTable(lines) {
  const cells = (line) => line.trim().replace(/^\|/, "").replace(/\|$/, "").split("|").map((cell) => cell.trim());
  const [header, , ...rows] = lines;
  const headHtml = cells(header).map((cell) => `<th>${inlineMarkdown(cell)}</th>`).join("");
  const bodyHtml = rows
    .map((row) => `<tr>${cells(row).map((cell) => `<td>${inlineMarkdown(cell)}</td>`).join("")}</tr>`)
    .join("");
  return `<div class="table-scroll"><table><thead><tr>${headHtml}</tr></thead><tbody>${bodyHtml}</tbody></table></div>`;
}

function isTableStart(lines, index) {
  return Boolean(
    lines[index] &&
      lines[index + 1] &&
      lines[index].includes("|") &&
      /^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)+\|?\s*$/.test(lines[index + 1])
  );
}

function renderMarkdownLite(markdown) {
  const lines = markdown.split(/\r?\n/);
  const html = [];
  let paragraph = [];
  let listType = "";
  let inCode = false;
  let code = [];

  const closeParagraph = () => {
    if (paragraph.length) {
      html.push(`<p>${inlineMarkdown(paragraph.join(" "))}</p>`);
      paragraph = [];
    }
  };

  const closeList = () => {
    if (listType) {
      html.push(`</${listType}>`);
      listType = "";
    }
  };

  for (let index = 0; index < lines.length; index += 1) {
    const line = lines[index];
    const trimmed = line.trim();

    if (trimmed.startsWith("```")) {
      closeParagraph();
      closeList();
      if (inCode) {
        html.push(`<pre><code>${escapeHtml(code.join("\n"))}</code></pre>`);
        code = [];
        inCode = false;
      } else {
        inCode = true;
      }
      continue;
    }

    if (inCode) {
      code.push(line);
      continue;
    }

    if (!trimmed) {
      closeParagraph();
      closeList();
      continue;
    }

    if (isTableStart(lines, index)) {
      closeParagraph();
      closeList();
      const tableLines = [lines[index], lines[index + 1]];
      index += 2;
      while (index < lines.length && lines[index].includes("|") && lines[index].trim()) {
        tableLines.push(lines[index]);
        index += 1;
      }
      index -= 1;
      html.push(renderTable(tableLines));
      continue;
    }

    if (/^<[^>]+>/.test(trimmed)) {
      closeParagraph();
      closeList();
      html.push(line);
      continue;
    }

    const heading = /^(#{1,6})\s+(.+)$/.exec(trimmed);
    if (heading) {
      closeParagraph();
      closeList();
      const level = heading[1].length;
      html.push(`<h${level}>${inlineMarkdown(heading[2])}</h${level}>`);
      continue;
    }

    if (trimmed.startsWith(">")) {
      closeParagraph();
      closeList();
      html.push(`<blockquote><p>${inlineMarkdown(trimmed.replace(/^>\s?/, ""))}</p></blockquote>`);
      continue;
    }

    const unordered = /^[-*]\s+(.+)$/.exec(trimmed);
    const ordered = /^\d+\.\s+(.+)$/.exec(trimmed);
    if (unordered || ordered) {
      closeParagraph();
      const nextType = unordered ? "ul" : "ol";
      if (listType !== nextType) {
        closeList();
        html.push(`<${nextType}>`);
        listType = nextType;
      }
      html.push(`<li>${inlineMarkdown((unordered || ordered)[1])}</li>`);
      continue;
    }

    closeList();
    paragraph.push(line);
  }

  closeParagraph();
  closeList();
  return html.join("\n");
}

function renderMarkdown(markdown, title = "") {
  const body = rewriteLegacyMarkdown(stripDuplicateTitleHeading(stripFrontMatter(markdown), title));
  if (window.marked && typeof window.marked.parse === "function") {
    window.marked.setOptions({ gfm: true, breaks: false });
    return window.marked.parse(body);
  }
  return renderMarkdownLite(body);
}

async function loadPostMarkdown(post, contentBySlug) {
  const embedded = cleanText(contentBySlug[post.slug]);
  if (embedded) {
    return embedded;
  }
  return fetchText(post.source);
}

function postMatches(post, query, tag) {
  const haystack = [
    post.title,
    post.date,
    post.topic,
    post.description,
    ...(post.tags || [])
  ].join(" ").toLowerCase();
  const matchesQuery = !query || haystack.includes(query.toLowerCase());
  const matchesTag = !tag || tag === "All" || (post.tags || []).includes(tag);
  return matchesQuery && matchesTag;
}

function createTagButton(tag, activeTag, onClick) {
  const button = document.createElement("button");
  button.className = "tag-button";
  button.type = "button";
  button.textContent = tag;
  button.setAttribute("aria-pressed", tag === activeTag ? "true" : "false");
  button.addEventListener("click", () => onClick(tag));
  return button;
}

function renderArchive(posts) {
  const list = document.getElementById("post-list");
  const filter = document.getElementById("tag-filter");
  const search = document.getElementById("post-search");
  if (!list || !filter || !search) {
    return;
  }

  let activeTag = "All";
  const tags = ["All", ...Array.from(new Set(posts.flatMap((post) => post.tags || []))).sort()];

  const paint = () => {
    const query = search.value.trim();
    const filtered = posts.filter((post) => postMatches(post, query, activeTag));
    list.replaceChildren();

    if (!filtered.length) {
      const empty = document.createElement("p");
      empty.className = "empty-state";
      empty.textContent = "No notes match this filter.";
      list.append(empty);
      return;
    }

    filtered.forEach((post) => {
      const article = document.createElement("article");
      article.className = "post-archive-item";
      article.innerHTML = `
        <div>
          <p class="section-kicker">${escapeHtml(post.topic || "Note")}</p>
          <h3>${escapeHtml(post.title)}</h3>
          <p>${escapeHtml(post.description || "")}</p>
          <div class="item-meta">
            <span>${escapeHtml(formatDate(post.date))}</span>
            ${(post.tags || []).map((tag) => `<span>${escapeHtml(tag)}</span>`).join("")}
          </div>
        </div>
        <a class="post-card-link" href="${escapeHtml(post.url)}" aria-label="Read ${escapeHtml(post.title)}"></a>
      `;
      list.append(article);
    });
  };

  const paintTags = () => {
    filter.replaceChildren();
    tags.forEach((tag) => filter.append(createTagButton(tag, activeTag, (nextTag) => {
      activeTag = nextTag;
      paintTags();
      paint();
    })));
  };

  search.addEventListener("input", paint);
  paintTags();
  paint();
}

function showMissingAssets(container) {
  container.querySelectorAll("img").forEach((image) => {
    image.addEventListener("error", () => {
      const note = document.createElement("p");
      note.className = "missing-asset";
      note.textContent = `Image unavailable: ${image.getAttribute("src")}`;
      image.replaceWith(note);
    }, { once: true });
  });
}

async function renderReader(posts) {
  const slug = getPostSlug();
  const post = posts.find((item) => item.slug === slug);
  const title = document.getElementById("post-title");
  const topic = document.getElementById("post-topic");
  const meta = document.getElementById("post-meta");
  const body = document.getElementById("post-body");

  if (!title || !topic || !meta || !body) {
    return;
  }

  if (!post) {
    document.title = "Note not found | Wang Jianghai";
    title.textContent = "Note not found";
    topic.textContent = "Archive";
    body.innerHTML = '<p class="empty-state">The requested note does not exist in the current post index.</p>';
    return;
  }

  document.title = `${post.title} | Wang Jianghai`;
  title.textContent = post.title;
  topic.textContent = post.topic || "Note";
  meta.innerHTML = `
    <span>${escapeHtml(formatDate(post.date))}</span>
    ${(post.tags || []).map((tag) => `<span>${escapeHtml(tag)}</span>`).join("")}
  `;

  try {
    const contentBySlug = await fetchJson(POSTS_CONTENT_URL);
    const markdown = await loadPostMarkdown(post, contentBySlug);
    await ensureMarked();
    body.innerHTML = renderMarkdown(markdown, post.title);
    showMissingAssets(body);
    await ensureMathJax();
    if (window.MathJax && window.MathJax.typesetPromise) {
      window.MathJax.typesetPromise([body]);
    }
  } catch (error) {
    body.innerHTML = `<p class="empty-state">${escapeHtml(error.message)}</p>`;
  }
}

document.addEventListener("DOMContentLoaded", async () => {
  try {
    const posts = await fetchJson(POSTS_INDEX_URL);
    const page = document.body.dataset.page;
    if (page === "posts-index") {
      renderArchive(posts);
    }
    if (page === "post-reader") {
      renderReader(posts);
    }
  } catch (error) {
    const target = document.getElementById("post-list") || document.getElementById("post-body");
    if (target) {
      target.innerHTML = `<p class="empty-state">${escapeHtml(error.message)}</p>`;
    }
  }
});
