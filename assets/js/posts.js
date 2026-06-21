const SITE_CONFIG_URL = "/assets/data/site.json";
const CONTENT_INDEX_URL = "/assets/data/content.json";
const POSTS_CONTENT_URL = "/assets/data/post-content.json";
const NEWS_CONTENT_URL = "/assets/data/news-content.json";
const MARKED_URL = "https://cdn.jsdelivr.net/npm/marked/marked.min.js";
const MATHJAX_URL = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js";

const cleanText = (value) => (value || "").toString().trim();
let currentLanguage = "en";

function setText(selector, value) {
  document.querySelectorAll(selector).forEach((node) => {
    node.textContent = cleanText(value);
  });
}

function setMetaContent(selector, value) {
  const node = document.querySelector(selector);
  if (node) {
    node.content = cleanText(value);
  }
}

function pageTitle(title, profileName) {
  return [cleanText(title), cleanText(profileName)].filter(Boolean).join(" | ");
}

function applySiteConfig(config, content, page) {
  const site = config.site || {};
  const profile = config.profile || {};
  const navigation = config.navigation || {};
  const postsPage = config.postsPage || {};
  const postPage = config.postPage || {};
  const newsPage = config.newsPage || {};

  currentLanguage = cleanText(site.language) || "en";
  document.documentElement.lang = currentLanguage;
  setText("[data-brand-mark]", site.brandMark);
  setText('[data-profile="name"]', profile.name);
  setText('[data-profile="location"]', profile.location);
  setText("[data-nav-news]", navigation.news);
  setText("[data-nav-research]", navigation.research);
  setText("[data-nav-writing]", navigation.writing);
  setMetaContent('meta[name="author"]', profile.name);

  const brand = document.querySelector("[data-brand-link]");
  if (brand) {
    brand.setAttribute("aria-label", `${cleanText(profile.name)} homepage`.trim());
  }

  document.querySelectorAll("[data-nav-news]").forEach((node) => {
    node.hidden = !Array.isArray(content.news) || !content.news.length;
  });

  if (page === "posts-index") {
    document.title = pageTitle(postsPage.title, profile.name);
    setMetaContent('meta[name="description"]', postsPage.description);
    ["eyebrow", "heading", "introduction", "searchLabel", "archiveKicker", "archiveTitle"].forEach((key) => {
      setText(`[data-posts-page="${key}"]`, postsPage[key]);
    });

    const toolbar = document.querySelector("[data-post-toolbar]");
    if (toolbar) {
      toolbar.setAttribute("aria-label", cleanText(postsPage.toolbarLabel));
    }
    const search = document.getElementById("post-search");
    if (search) {
      search.placeholder = cleanText(postsPage.searchPlaceholder);
    }
    const filter = document.getElementById("tag-filter");
    if (filter) {
      filter.setAttribute("aria-label", cleanText(postsPage.tagFilterLabel));
    }
  }

  if (page === "post-reader") {
    document.title = pageTitle(postPage.title, profile.name);
    setMetaContent('meta[name="description"]', postPage.description);
    setText('[data-post-page="backLabel"]', postPage.backLabel);
    setText("#post-topic", postPage.defaultTopic);
    setText("#post-title", postPage.loadingTitle);
  }

  if (page === "news-reader") {
    document.title = pageTitle(newsPage.title, profile.name);
    setMetaContent('meta[name="description"]', newsPage.description);
    setText('[data-news-page="backLabel"]', newsPage.backLabel);
    setText("#news-topic", newsPage.defaultTopic);
    setText("#news-title", newsPage.loadingTitle);
  }
}

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
  return new Intl.DateTimeFormat(currentLanguage, {
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
  const matchesTag = !tag || (post.tags || []).includes(tag);
  return matchesQuery && matchesTag;
}

function createTagButton(value, label, activeTag, onClick) {
  const button = document.createElement("button");
  button.className = "tag-button";
  button.type = "button";
  button.textContent = label;
  button.setAttribute("aria-pressed", value === activeTag ? "true" : "false");
  button.addEventListener("click", () => onClick(value));
  return button;
}

function renderArchive(posts, config) {
  const postsPage = config.postsPage || {};
  const postPage = config.postPage || {};
  const list = document.getElementById("post-list");
  const filter = document.getElementById("tag-filter");
  const search = document.getElementById("post-search");
  if (!list || !filter || !search) {
    return;
  }

  let activeTag = "";
  const tags = ["", ...Array.from(new Set(posts.flatMap((post) => post.tags || []))).sort()];

  const paint = () => {
    const query = search.value.trim();
    const filtered = posts.filter((post) => postMatches(post, query, activeTag));
    list.replaceChildren();

    if (!filtered.length) {
      const empty = document.createElement("p");
      empty.className = "empty-state";
      empty.textContent = cleanText(postsPage.emptyMessage);
      list.append(empty);
      return;
    }

    filtered.forEach((post) => {
      const article = document.createElement("article");
      article.className = "post-archive-item";
      article.innerHTML = `
        <div>
          <p class="section-kicker">${escapeHtml(post.topic || postPage.defaultTopic)}</p>
          <h3>${escapeHtml(post.title)}</h3>
          <p>${escapeHtml(post.description || "")}</p>
          <div class="item-meta">
            <span>${escapeHtml(formatDate(post.date))}</span>
            ${(post.tags || []).map((tag) => `<span>${escapeHtml(tag)}</span>`).join("")}
          </div>
        </div>
        <a class="post-card-link" href="${escapeHtml(post.url)}" aria-label="${escapeHtml(postPage.readLabel)} ${escapeHtml(post.title)}"></a>
      `;
      list.append(article);
    });
  };

  const paintTags = () => {
    filter.replaceChildren();
    tags.forEach((tag) => filter.append(createTagButton(
      tag,
      tag || cleanText(postsPage.allTagsLabel),
      activeTag,
      (nextTag) => {
      activeTag = nextTag;
      paintTags();
      paint();
      }
    )));
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

async function renderReader(posts, config) {
  const profile = config.profile || {};
  const postPage = config.postPage || {};
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
    document.title = pageTitle(postPage.notFoundTitle, profile.name);
    title.textContent = cleanText(postPage.notFoundTitle);
    topic.textContent = cleanText(postPage.notFoundTopic);
    body.innerHTML = `<p class="empty-state">${escapeHtml(postPage.notFoundMessage)}</p>`;
    return;
  }

  document.title = pageTitle(post.title, profile.name);
  title.textContent = post.title;
  topic.textContent = post.topic || cleanText(postPage.defaultTopic);
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

async function renderNewsReader(news, config) {
  const profile = config.profile || {};
  const newsPage = config.newsPage || {};
  const slug = getPostSlug();
  const item = news.find((entry) => entry.slug === slug);
  const title = document.getElementById("news-title");
  const topic = document.getElementById("news-topic");
  const meta = document.getElementById("news-meta");
  const image = document.getElementById("news-image");
  const summary = document.getElementById("news-summary");
  const body = document.getElementById("news-body");

  if (!title || !topic || !meta || !image || !summary || !body) {
    return;
  }

  if (!item) {
    document.title = pageTitle(newsPage.notFoundTitle, profile.name);
    title.textContent = cleanText(newsPage.notFoundTitle);
    topic.textContent = cleanText(newsPage.notFoundTopic);
    body.innerHTML = `<p class="empty-state">${escapeHtml(newsPage.notFoundMessage)}</p>`;
    return;
  }

  document.title = pageTitle(item.title, profile.name);
  title.textContent = item.title;
  topic.textContent = cleanText(newsPage.defaultTopic);
  meta.innerHTML = `<span>${escapeHtml(formatDate(item.date))}</span>`;
  image.src = cleanText(item.image);
  image.alt = cleanText(item.imageAlt);
  image.hidden = !image.src;
  summary.textContent = cleanText(item.description);
  summary.hidden = !summary.textContent;
  if (!image.hidden && image.parentElement) {
    showMissingAssets(image.parentElement);
  }

  try {
    const contentBySlug = await fetchJson(NEWS_CONTENT_URL);
    const markdown = await loadPostMarkdown(item, contentBySlug);
    await ensureMarked();
    body.innerHTML = renderMarkdown(markdown, item.title);
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
    const [config, content] = await Promise.all([
      fetchJson(SITE_CONFIG_URL),
      fetchJson(CONTENT_INDEX_URL)
    ]);
    const posts = Array.isArray(content.posts) ? content.posts : [];
    const news = Array.isArray(content.news) ? content.news : [];
    const page = document.body.dataset.page;
    applySiteConfig(config, content, page);
    if (page === "posts-index") {
      renderArchive(posts, config);
    }
    if (page === "post-reader") {
      renderReader(posts, config);
    }
    if (page === "news-reader") {
      renderNewsReader(news, config);
    }
  } catch (error) {
    const target = document.getElementById("post-list") || document.getElementById("post-body") || document.getElementById("news-body");
    if (target) {
      target.innerHTML = `<p class="empty-state">${escapeHtml(error.message)}</p>`;
    }
  }
});
