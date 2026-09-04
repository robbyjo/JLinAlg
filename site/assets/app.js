const THEME_KEY = "jlinalg-theme";

function applyTheme(value) {
  const theme = value === "light" || value === "dark" ? value : "system";
  if (theme === "system") document.documentElement.removeAttribute("data-theme");
  else document.documentElement.dataset.theme = theme;
  document.querySelectorAll("[data-theme-select]").forEach(select => {
    select.value = theme;
  });
}

try { applyTheme(localStorage.getItem(THEME_KEY) || "system"); } catch { applyTheme("system"); }

class SiteHeader extends HTMLElement {
  connectedCallback() {
    const root = this.getAttribute("root") || "./";
    this.innerHTML = `
      <a class="skip-link" href="#main">Skip to content</a>
      <header class="site-header">
        <div class="container nav-wrap">
          <a class="brand" href="${root}index.html" aria-label="JLinAlg home">
            <span class="brand-mark" aria-hidden="true"><i></i><i></i><i></i><i></i><i></i><i></i></span>
            <span>JLinAlg</span>
          </a>
          <button class="nav-toggle" type="button" aria-label="Open navigation" aria-expanded="false">☰</button>
          <nav class="nav-links" aria-label="Primary navigation">
            <a href="${root}index.html#features">Features</a>
            <a href="${root}vignettes/index.html">Vignettes</a>
            <a href="${root}verification.html">Verification</a>
            <a href="${root}index.html#architecture">Architecture</a>
          </nav>
          <div class="nav-actions">
            <label>
              <span class="skip-link">Color scheme</span>
              <select class="theme-select" data-theme-select aria-label="Color scheme">
                <option value="system">System</option><option value="light">Light</option><option value="dark">Dark</option>
              </select>
            </label>
            <a class="button secondary" href="https://github.com/robbyjo/JLinAlg">GitHub ↗</a>
          </div>
        </div>
      </header>`;
    const toggle = this.querySelector(".nav-toggle");
    const nav = this.querySelector(".nav-links");
    toggle.addEventListener("click", () => {
      const open = nav.classList.toggle("open");
      toggle.setAttribute("aria-expanded", String(open));
      toggle.textContent = open ? "×" : "☰";
    });
    this.querySelector("[data-theme-select]").addEventListener("change", event => {
      const value = event.target.value;
      try { localStorage.setItem(THEME_KEY, value); } catch { /* private mode */ }
      applyTheme(value);
    });
    applyTheme((() => { try { return localStorage.getItem(THEME_KEY) || "system"; } catch { return "system"; } })());
  }
}

class SiteFooter extends HTMLElement {
  connectedCallback() {
    const root = this.getAttribute("root") || "./";
    this.innerHTML = `
      <footer class="site-footer">
        <div class="container">
          <div class="footer-grid">
            <div><a class="brand" href="${root}index.html"><span class="brand-mark" aria-hidden="true"><i></i><i></i><i></i><i></i><i></i><i></i></span><span>JLinAlg</span></a><p>High-performance Java statistical models for genetic association, causal inference, fine mapping, and quantitative research.</p></div>
            <div class="footer-links"><strong>Learn</strong><a href="${root}vignettes/index.html">All vignettes</a><a href="${root}verification.html">Accuracy and performance</a><a href="https://github.com/robbyjo/JLinAlg/blob/main/docs/numerical-contract.md">Numerical contract</a></div>
            <div class="footer-links"><strong>Project</strong><a href="https://github.com/robbyjo/JLinAlg">Source code</a><a href="https://github.com/robbyjo/JLinAlg/issues">Issues</a><a href="https://github.com/robbyjo/JDistlib">JDistlib</a></div>
          </div>
          <div class="footer-bottom"><span>© <span data-year></span> JLinAlg contributors · GPL-2.0-or-later</span><span>JLinAlg v0.2.0 · Java 17+ · JDistlib 0.10.1 · FP64</span></div>
        </div>
      </footer>`;
    this.querySelector("[data-year]").textContent = new Date().getFullYear();
  }
}

customElements.define("site-header", SiteHeader);
customElements.define("site-footer", SiteFooter);

document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll(".code-window").forEach(window => {
    const bar = window.querySelector(".code-bar");
    const code = window.querySelector("code");
    if (!bar || !code || bar.querySelector(".copy-code")) return;
    const button = document.createElement("button");
    button.className = "copy-code";
    button.type = "button";
    button.textContent = "Copy";
    button.addEventListener("click", async () => {
      try {
        await navigator.clipboard.writeText(code.innerText);
        button.textContent = "Copied";
        setTimeout(() => { button.textContent = "Copy"; }, 1500);
      } catch { button.textContent = "Select code"; }
    });
    bar.append(button);
  });

  const filters = document.querySelectorAll("[data-filter]");
  const cards = document.querySelectorAll("[data-category]");
  filters.forEach(button => button.addEventListener("click", () => {
    filters.forEach(item => item.setAttribute("aria-pressed", "false"));
    button.setAttribute("aria-pressed", "true");
    const filter = button.dataset.filter;
    cards.forEach(card => { card.hidden = filter !== "all" && card.dataset.category !== filter; });
  }));

  const reduced = matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (!reduced && "IntersectionObserver" in window) {
    const observer = new IntersectionObserver(entries => entries.forEach(entry => {
      if (entry.isIntersecting) { entry.target.classList.add("visible"); observer.unobserve(entry.target); }
    }), { threshold: .08 });
    document.querySelectorAll(".reveal").forEach(item => observer.observe(item));
  } else document.querySelectorAll(".reveal").forEach(item => item.classList.add("visible"));

  const tocLinks = [...document.querySelectorAll(".toc a")];
  if (tocLinks.length && "IntersectionObserver" in window) {
    const sections = tocLinks.map(link => document.querySelector(link.getAttribute("href"))).filter(Boolean);
    const tocObserver = new IntersectionObserver(entries => entries.forEach(entry => {
      if (entry.isIntersecting) {
        tocLinks.forEach(link => link.classList.toggle("active", link.getAttribute("href") === `#${entry.target.id}`));
      }
    }), { rootMargin: "-20% 0px -70% 0px" });
    sections.forEach(section => tocObserver.observe(section));
  }
});
