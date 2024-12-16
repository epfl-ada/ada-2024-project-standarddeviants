document.addEventListener("DOMContentLoaded", () => {
    const navbar = document.getElementById("custom-navbar");
    if (navbar) {
      const navbarHeight = 1. * navbar.offsetHeight;
      document.documentElement.style.setProperty('--navbar-height', `${navbarHeight}px`);
    }
  });
