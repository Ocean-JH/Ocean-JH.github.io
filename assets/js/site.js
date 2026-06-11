const defaultSiteData = {
  profile: {
    name: "Wang Jianghai",
    nativeName: "王江海",
    role: "PhD student in Computational Materials Science",
    affiliation: "Nanyang Technological University",
    location: "Singapore",
    email: "jianghai001@e.ntu.edu.sg",
    github: "Ocean-JH",
    cvUrl: ""
  },
  researchInterests: [
    {
      title: "Crystal structure prediction",
      description: "Exploring metastable materials through global optimization, structural prototypes, and first-principles validation."
    },
    {
      title: "First-principles workflows",
      description: "Connecting electronic structure, phonons, elasticity, dielectric response, and transport calculations for materials discovery."
    },
    {
      title: "AI for materials",
      description: "Using machine learning potentials and data-driven models to accelerate configuration-space exploration."
    }
  ],
  publications: [],
  writing: []
};

const text = (value) => (value || "").toString().trim();

function setText(selector, value) {
  document.querySelectorAll(selector).forEach((node) => {
    node.textContent = value;
  });
}

function setHref(selector, href) {
  document.querySelectorAll(selector).forEach((node) => {
    node.href = href;
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
  return {
    ...defaultSiteData,
    ...data,
    profile: {
      ...defaultSiteData.profile,
      ...(data.profile || {})
    },
    researchInterests: Array.isArray(data.researchInterests)
      ? data.researchInterests
      : defaultSiteData.researchInterests,
    publications: Array.isArray(data.publications) ? data.publications : [],
    writing: Array.isArray(data.writing) ? data.writing : []
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

async function loadPostsData() {
  try {
    const response = await fetch("/assets/data/posts.json", { cache: "no-store" });
    if (!response.ok) {
      throw new Error(`Unable to load posts data: ${response.status}`);
    }
    const posts = await response.json();
    return Array.isArray(posts) ? posts : [];
  } catch (error) {
    console.warn(error.message);
    return [];
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

function updateProfile(profile) {
  setText('[data-profile="name"]', text(profile.name));
  setText('[data-profile="nativeName"]', text(profile.nativeName));
  setText('[data-profile="role"]', text(profile.role));
  setText('[data-profile="affiliation"]', text(profile.affiliation));
  setText('[data-profile="location"]', text(profile.location));

  const email = text(profile.email);
  if (email) {
    setHref("[data-email-link]", `mailto:${email}`);
    setText("[data-email-text]", email);
  }

  const github = text(profile.github);
  if (github) {
    setHref("[data-github-link]", `https://github.com/${github}`);
    setText("[data-github-text]", `github.com/${github}`);
  }

  document.querySelectorAll("[data-cv-link]").forEach((node) => {
    const cvUrl = text(profile.cvUrl);
    node.hidden = !cvUrl;
    if (cvUrl) {
      node.href = cvUrl;
    }
  });
}

function renderResearch(items) {
  const section = document.querySelector('[data-section="research"]');
  const list = document.getElementById("research-list");
  if (!section || !list) {
    return;
  }

  if (!items.length) {
    section.hidden = true;
    return;
  }

  list.replaceChildren();
  items.forEach((item) => {
    const card = createElement("article", "interest-card");
    card.append(
      createElement("h3", "", text(item.title)),
      createElement("p", "", text(item.description))
    );
    list.append(card);
  });
  section.hidden = false;
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
    const article = createElement("article", "writing-item");
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

    const link = createSafeLink("Read", text(item.url));
    if (link) {
      const links = createElement("div", "item-links");
      links.append(link);
      article.append(links);
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

document.addEventListener("DOMContentLoaded", async () => {
  const [siteData, posts] = await Promise.all([loadSiteData(), loadPostsData()]);
  updateProfile(siteData.profile);
  renderResearch(siteData.researchInterests);
  renderPublications(siteData.publications);
  const writingItems = siteData.writing.length
    ? siteData.writing
    : posts.slice(0, 6).map(postToWritingItem);
  renderWriting(writingItems, posts.length > 0);
});
