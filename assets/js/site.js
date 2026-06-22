const defaultSiteData = {
  site: {
    language: "en",
    url: "",
    title: "",
    description: "",
    socialDescription: "",
    brandMark: ""
  },
  navigation: {
    news: "",
    research: "",
    publications: "",
    writing: ""
  },
  profile: {
    name: "",
    nativeName: "",
    role: "",
    affiliation: "",
    location: "",
    email: "",
    github: "",
    linkedin: "",
    orcid: "",
    googleScholar: "",
    cvUrl: ""
  },
  home: {
    hero: {
      eyebrow: "",
      description: "",
      socialAriaLabel: "",
      researchMarkLabel: "",
      researchMarkCaptions: []
    },
    summary: {
      kicker: "",
      title: "",
      description: ""
    },
    news: { kicker: "", title: "" },
    research: { kicker: "", title: "" },
    publications: { kicker: "", title: "" },
    writing: {
      kicker: "",
      title: "",
      archiveLabel: "",
      archiveUrl: "/posts/",
      homepageLimit: 6
    }
  },
  postsPage: {},
  postPage: {},
  resources: []
};

const defaultContentData = {
  news: [],
  research: [],
  publications: [],
  posts: []
};

const text = (value) => (value || "").toString().trim();

function setText(selector, value) {
  document.querySelectorAll(selector).forEach((node) => {
    node.textContent = value;
  });
}

function setMetaContent(selector, value) {
  const node = document.querySelector(selector);
  if (node) {
    node.content = text(value);
  }
}

function updateOptionalLink(selector, href) {
  document.querySelectorAll(selector).forEach((node) => {
    node.hidden = !href;
    if (href) {
      node.href = href;
    }
  });
}

function createElement(tag, className, content) {
  const element = document.createElement(tag);
  if (className) {
    element.className = className;
  }
  if (content) {
    element.textContent = content;
  }
  return element;
}

function createSafeLink(label, href) {
  if (!href) {
    return null;
  }
  const link = createElement("a", "", label);
  link.href = href;
  link.rel = "noopener";
  return link;
}

function doiUrl(doi) {
  const value = text(doi);
  if (!value) {
    return "";
  }
  return value.startsWith("http") ? value : `https://doi.org/${value}`;
}

function mergeSiteData(data) {
  const home = data.home || {};
  return {
    ...defaultSiteData,
    ...data,
    site: {
      ...defaultSiteData.site,
      ...(data.site || {})
    },
    navigation: {
      ...defaultSiteData.navigation,
      ...(data.navigation || {})
    },
    profile: {
      ...defaultSiteData.profile,
      ...(data.profile || {})
    },
    home: {
      ...defaultSiteData.home,
      ...home,
      hero: {
        ...defaultSiteData.home.hero,
        ...(home.hero || {}),
        researchMarkCaptions: Array.isArray(home.hero && home.hero.researchMarkCaptions)
          ? home.hero.researchMarkCaptions
          : []
      },
      summary: {
        ...defaultSiteData.home.summary,
        ...(home.summary || {})
      },
      news: {
        ...defaultSiteData.home.news,
        ...(home.news || {})
      },
      research: {
        ...defaultSiteData.home.research,
        ...(home.research || {})
      },
      publications: {
        ...defaultSiteData.home.publications,
        ...(home.publications || {})
      },
      writing: {
        ...defaultSiteData.home.writing,
        ...(home.writing || {})
      },
      resources: {
        kicker: "",
        title: "",
        ...(home.resources || {})
      }
    },
    postsPage: {
      ...(data.postsPage || {})
    },
    postPage: {
      ...(data.postPage || {})
    },
    resources: Array.isArray(data.resources) ? data.resources : []
  };
}

async function loadSiteData() {
  try {
    const response = await fetch("/assets/data/site.json", { cache: "no-store" });
    if (!response.ok) {
      throw new Error(`Unable to load site data: ${response.status}`);
    }
    return mergeSiteData(await response.json());
  } catch (error) {
    console.warn(error.message);
    return defaultSiteData;
  }
}

async function loadContentData() {
  try {
    const response = await fetch("/assets/data/content.json", { cache: "no-store" });
    if (!response.ok) {
      throw new Error(`Unable to load content data: ${response.status}`);
    }
    const data = await response.json();
    return {
      news: Array.isArray(data.news) ? data.news : [],
      research: Array.isArray(data.research) ? data.research : [],
      publications: Array.isArray(data.publications) ? data.publications : [],
      posts: Array.isArray(data.posts) ? data.posts : []
    };
  } catch (error) {
    console.warn(error.message);
    return defaultContentData;
  }
}

function postToWritingItem(post) {
  return {
    title: post.title,
    description: post.description,
    date: post.date,
    topic: post.topic,
    url: post.url
  };
}

function formatDate(value) {
  const parts = text(value).split("-").map(Number);
  if (parts.length !== 3 || parts.some((part) => !part)) {
    return text(value);
  }
  return new Intl.DateTimeFormat(document.documentElement.lang || "en", {
    year: "numeric",
    month: "short",
    day: "numeric"
  }).format(new Date(Date.UTC(parts[0], parts[1] - 1, parts[2])));
}

function updateSiteMetadata(site, profile) {
  const language = text(site.language) || "en";
  document.documentElement.lang = language;

  if (text(site.title)) {
    document.title = text(site.title);
  }
  setMetaContent('meta[name="description"]', site.description);
  setMetaContent('meta[name="author"]', profile.name);
  setMetaContent('meta[property="og:title"]', site.title);
  setMetaContent('meta[property="og:description"]', site.socialDescription || site.description);
  setMetaContent('meta[property="og:url"]', site.url);

  setText("[data-brand-mark]", text(site.brandMark));
  const brand = document.querySelector("[data-brand-link]");
  if (brand) {
    brand.setAttribute("aria-label", `${text(profile.name)} homepage`.trim());
  }
}

function updateNavigation(navigation) {
  setText("[data-nav-news]", text(navigation.news));
  setText("[data-nav-research]", text(navigation.research));
  setText("[data-nav-publications]", text(navigation.publications));
  setText("[data-nav-writing]", text(navigation.writing));
}

function updateHomeCopy(home) {
  const values = {
    heroEyebrow: home.hero.eyebrow,
    heroDescription: home.hero.description,
    summaryKicker: home.summary.kicker,
    summaryTitle: home.summary.title,
    summaryDescription: home.summary.description,
    newsKicker: home.news.kicker,
    newsTitle: home.news.title,
    researchKicker: home.research.kicker,
    researchTitle: home.research.title,
    publicationsKicker: home.publications.kicker,
    publicationsTitle: home.publications.title,
    writingKicker: home.writing.kicker,
    writingTitle: home.writing.title,
    resourcesKicker: home.resources.kicker,
    resourcesTitle: home.resources.title
  };

  Object.entries(values).forEach(([key, value]) => {
    setText(`[data-home="${key}"]`, text(value));
  });

  const socials = document.querySelector("[data-hero-socials]");
  if (socials) {
    socials.setAttribute("aria-label", text(home.hero.socialAriaLabel));
  }

  const researchMark = document.querySelector("[data-research-mark]");
  if (researchMark) {
    researchMark.setAttribute("aria-label", text(home.hero.researchMarkLabel));
  }

  const captions = document.getElementById("research-mark-caption");
  if (captions) {
    captions.replaceChildren();
    home.hero.researchMarkCaptions.map(text).filter(Boolean).forEach((caption) => {
      captions.append(createElement("span", "", caption));
    });
  }

  const archiveLink = document.querySelector("[data-writing-archive-link]");
  if (archiveLink) {
    archiveLink.textContent = text(home.writing.archiveLabel);
    archiveLink.href = text(home.writing.archiveUrl) || "/posts/";
  }
}

function updateProfile(profile) {
  setText('[data-profile="name"]', text(profile.name));
  setText('[data-profile="nativeName"]', text(profile.nativeName));
  setText('[data-profile="role"]', text(profile.role));
  setText('[data-profile="affiliation"]', text(profile.affiliation));
  setText('[data-profile="location"]', text(profile.location));

  const email = text(profile.email);
  updateOptionalLink("[data-email-link]", email ? `mailto:${email}` : "");

  const github = text(profile.github);
  const githubUrl = github.startsWith("http") ? github : github ? `https://github.com/${github}` : "";
  updateOptionalLink("[data-github-link]", githubUrl);
  updateOptionalLink("[data-linkedin-link]", text(profile.linkedin));
  updateOptionalLink("[data-orcid-link]", text(profile.orcid));
  updateOptionalLink("[data-google-scholar-link]", text(profile.googleScholar));

  document.querySelectorAll("[data-cv-link]").forEach((node) => {
    const cvUrl = text(profile.cvUrl);
    node.hidden = !cvUrl;
    if (cvUrl) {
      node.href = cvUrl;
    }
  });
}

function renderNews(items) {
  const section = document.querySelector('[data-section="news"]');
  const nav = document.querySelector("[data-nav-news]");
  const list = document.getElementById("news-list");
  if (!section || !list) {
    return;
  }

  if (!items.length) {
    section.hidden = true;
    if (nav) {
      nav.hidden = true;
    }
    return;
  }

  list.replaceChildren();
  items.forEach((item) => {
    const url = text(item.url);
    const entry = createElement(url ? "a" : "article", url ? "news-item news-item-link" : "news-item");
    if (url) {
      entry.href = url;
      entry.rel = "noopener";
    }

    const date = createElement("time", "news-date", formatDate(item.date));
    if (text(item.date)) {
      date.dateTime = text(item.date);
    }

    const content = createElement("div", "news-content");
    content.append(createElement("h3", "", text(item.title)));
    entry.append(date, content);

    if (text(item.image)) {
      entry.classList.add("news-item-has-image");
      const image = createElement("img", "news-image");
      image.src = text(item.image);
      image.alt = text(item.imageAlt);
      image.loading = "lazy";
      image.decoding = "async";
      entry.append(image);
    }
    list.append(entry);
  });

  section.hidden = false;
  if (nav) {
    nav.hidden = false;
  }
}

function renderResearch(items) {
  const section = document.querySelector('[data-section="research"]');
  const nav = document.querySelector("[data-nav-research]");
  const list = document.getElementById("research-list");
  if (!section || !list) {
    return;
  }

  if (!items.length) {
    section.hidden = true;
    if (nav) {
      nav.hidden = true;
    }
    return;
  }

  list.replaceChildren();
  items.forEach((item) => {
    const url = text(item.url);
    const card = createElement(url ? "a" : "article", url ? "interest-card interest-card-link" : "interest-card");
    if (url) {
      card.href = url;
      card.rel = "noopener";
    }
    card.append(
      createElement("h3", "", text(item.title)),
      createElement("p", "", text(item.description))
    );
    list.append(card);
  });
  section.hidden = false;
  if (nav) {
    nav.hidden = false;
  }
}

function renderPublications(items) {
  const section = document.querySelector('[data-section="publications"]');
  const nav = document.querySelector("[data-nav-publications]");
  const list = document.getElementById("publication-list");
  if (!section || !list) {
    return;
  }

  const selected = items.some((item) => item.selected)
    ? items.filter((item) => item.selected)
    : items;

  if (!selected.length) {
    section.hidden = true;
    if (nav) {
      nav.hidden = true;
    }
    return;
  }

  list.replaceChildren();
  selected.forEach((item) => {
    const article = createElement("article", "publication-item");
    article.append(createElement("h3", "", text(item.title)));

    const authors = text(item.authors);
    if (authors) {
      article.append(createElement("p", "", authors));
    }

    const meta = createElement("div", "item-meta");
    [item.venue, item.year].map(text).filter(Boolean).forEach((value) => {
      meta.append(createElement("span", "", value));
    });
    if (meta.children.length) {
      article.append(meta);
    }

    const links = createElement("div", "item-links");
    [
      createSafeLink("DOI", doiUrl(item.doi)),
      createSafeLink("Article", text(item.url)),
      createSafeLink("PDF", text(item.pdf))
    ].filter(Boolean).forEach((link) => links.append(link));
    if (links.children.length) {
      article.append(links);
    }

    list.append(article);
  });

  section.hidden = false;
  if (nav) {
    nav.hidden = false;
  }
}

function renderWriting(items, showArchiveLink = false) {
  const section = document.querySelector('[data-section="writing"]');
  const nav = document.querySelector("[data-nav-writing]");
  const list = document.getElementById("writing-list");
  const actions = document.querySelector("[data-writing-actions]");
  if (!section || !list) {
    return;
  }

  if (!items.length) {
    section.hidden = true;
    if (actions) {
      actions.hidden = true;
    }
    if (nav) {
      nav.hidden = true;
    }
    return;
  }

  list.replaceChildren();
  items.forEach((item) => {
    const url = text(item.url);
    const article = createElement(url ? "a" : "article", url ? "writing-item writing-card-link" : "writing-item");
    if (url) {
      article.href = url;
      article.setAttribute("aria-label", `Read ${text(item.title)}`);
    }
    article.append(createElement("h3", "", text(item.title)));

    const description = text(item.description);
    if (description) {
      article.append(createElement("p", "", description));
    }

    const meta = createElement("div", "item-meta");
    [item.date, item.topic].map(text).filter(Boolean).forEach((value) => {
      meta.append(createElement("span", "", value));
    });
    if (meta.children.length) {
      article.append(meta);
    }

    list.append(article);
  });

  section.hidden = false;
  if (actions) {
    actions.hidden = !showArchiveLink;
  }
  if (nav) {
    nav.hidden = false;
  }
}

function renderResources(items) {
  const section = document.querySelector('[data-section="resources"]');
  const list = document.getElementById("resource-list");
  if (!section || !list) {
    return;
  }

  if (!items.length) {
    section.hidden = true;
    return;
  }

  list.replaceChildren();
  items.forEach((item) => {
    const url = text(item.url);
    const card = createElement(url ? "a" : "article", "resource-card");
    if (url) {
      card.href = url;
      card.rel = "noopener";
    }
    card.append(
      createElement("h3", "", text(item.title)),
      createElement("p", "", text(item.description))
    );
    list.append(card);
  });
  section.hidden = false;
}

document.addEventListener("DOMContentLoaded", async () => {
  const [siteData, contentData] = await Promise.all([loadSiteData(), loadContentData()]);
  updateSiteMetadata(siteData.site, siteData.profile);
  updateNavigation(siteData.navigation);
  updateProfile(siteData.profile);
  updateHomeCopy(siteData.home);
  renderNews(contentData.news);
  renderResearch(contentData.research);
  renderPublications(contentData.publications);
  const configuredLimit = Number(siteData.home.writing.homepageLimit);
  const homepageLimit = Number.isFinite(configuredLimit) && configuredLimit > 0 ? configuredLimit : 6;
  const writingItems = contentData.posts.slice(0, homepageLimit).map(postToWritingItem);
  renderWriting(writingItems, contentData.posts.length > 0 && Boolean(text(siteData.home.writing.archiveUrl)));
  renderResources(siteData.resources);
});
