// Copy-to-clipboard for the hero install command (and any [data-copy-target]).
(function () {
  function fallbackCopy(text) {
    var ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.style.position = "fixed";
    ta.style.top = "-9999px";
    document.body.appendChild(ta);
    ta.select();
    try { document.execCommand("copy"); } catch (e) { /* ignore */ }
    document.body.removeChild(ta);
  }

  // Icon-only feedback: swap to a checkmark briefly via CSS class (never touch
  // textContent here, or it would wipe out the button's SVG icon).
  function flash(btn) {
    clearTimeout(btn._copyTimer);
    btn.classList.add("is-copied");
    btn._copyTimer = setTimeout(function () {
      btn.classList.remove("is-copied");
    }, 1400);
  }

  // Straight-quote the displayed command from its verbatim data-cmd (pandoc's
  // smart-quotes mangle the body text, but leave raw-HTML attributes alone).
  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll("[data-cmd]").forEach(function (el) {
      el.textContent = el.getAttribute("data-cmd");
    });
  });

  document.addEventListener("click", function (e) {
    var btn = e.target.closest("[data-copy-target]");
    if (!btn) return;
    var el = document.getElementById(btn.getAttribute("data-copy-target"));
    if (!el) return;
    var text = el.getAttribute("data-cmd") || el.textContent;

    if (navigator.clipboard && navigator.clipboard.writeText) {
      // Guard against a permission prompt that never resolves (rare, but would
      // otherwise leave the button silently unresponsive): fall back after 1s.
      var settled = false;
      var give_up = setTimeout(function () {
        if (settled) return;
        settled = true;
        fallbackCopy(text);
        flash(btn);
      }, 1000);
      navigator.clipboard.writeText(text).then(
        function () { if (settled) return; settled = true; clearTimeout(give_up); flash(btn); },
        function () { if (settled) return; settled = true; clearTimeout(give_up); fallbackCopy(text); flash(btn); }
      );
    } else {
      fallbackCopy(text);
      flash(btn);
    }
  });
})();
