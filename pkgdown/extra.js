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

  function flash(btn) {
    var prev = btn.getAttribute("data-label") || btn.textContent;
    btn.setAttribute("data-label", prev);
    btn.textContent = "Copied";
    btn.classList.add("is-copied");
    setTimeout(function () {
      btn.textContent = prev;
      btn.classList.remove("is-copied");
    }, 1600);
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
      navigator.clipboard.writeText(text).then(
        function () { flash(btn); },
        function () { fallbackCopy(text); flash(btn); }
      );
    } else {
      fallbackCopy(text);
      flash(btn);
    }
  });
})();
