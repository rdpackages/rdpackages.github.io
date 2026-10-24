(function () {
  var menuToggle = document.getElementById("site-nav-toggle");
  var navigation = document.getElementById("site-nav");
  var mobileLayout = window.matchMedia("(max-width: 840px)");

  if (!menuToggle || !navigation) {
    return;
  }

  function closeMenu(restoreFocus) {
    menuToggle.checked = false;
    if (restoreFocus && mobileLayout.matches) {
      menuToggle.focus();
    }
  }

  navigation.querySelectorAll("a").forEach(function (link) {
    link.addEventListener("click", function () {
      closeMenu();
      var destination = new URL(link.href);
      if (destination.pathname === window.location.pathname && destination.hash) {
        var target = document.getElementById(decodeURIComponent(destination.hash.slice(1)));
        if (target) {
          if (!target.hasAttribute("tabindex")) {
            target.setAttribute("tabindex", "-1");
            target.addEventListener("blur", function () {
              target.removeAttribute("tabindex");
            }, { once: true });
          }
          target.focus({ preventScroll: true });
        }
      }
    });
  });

  document.addEventListener("keydown", function (event) {
    if (event.key === "Escape" && menuToggle.checked) {
      closeMenu(true);
    }
  });

  document.addEventListener("click", function (event) {
    if (!menuToggle.parentElement.contains(event.target)) {
      closeMenu();
    }
  });

  menuToggle.parentElement.addEventListener("focusout", function (event) {
    if (event.relatedTarget && !menuToggle.parentElement.contains(event.relatedTarget)) {
      closeMenu();
    }
  });

  mobileLayout.addEventListener("change", function () {
    closeMenu(navigation.contains(document.activeElement));
  });

  window.addEventListener("pagehide", function () { closeMenu(); });
  window.addEventListener("pageshow", function () { closeMenu(); });
})();
