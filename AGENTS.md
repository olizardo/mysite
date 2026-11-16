# Agent Instructions

This is the repository for my academic website, built using Quarto.
When interacting with this codebase, please adhere to the following rules:

1. **Quarto:** Remember that this site relies on Quarto. Any changes to navigation should be made in `_quarto.yml`. 
2. **Re-rendering:** Run `quarto render` or `quarto render <file/folder>` to build and test changes to the website structure or content.
3. **Data/Content Migrations:** Keep consistent formatting for migrated content (e.g. keeping YAML frontmatter clean and minimal, using the standard `.qmd` extensions).
4. **CV Generation:** The CV is a single, self-contained Quarto document. Any updates, additions, or modifications to the CV should be made directly in `cv.qmd` using standard Markdown formatting.
5. **Themes/CSS:** For aesthetic or visual adjustments, utilize the custom stylesheets (`custom.scss`, `custom-dark.scss`, and `styles.css`).

## Content conventions (2026-11-03)

- **cv.qmd:** Publications are organized as `##` sections (Books, Peer-Reviewed Articles, Book Chapters, Comments & Short Pieces, Encyclopedia Entries, etc.) with `###` period subheadings; each entry is one paragraph containing a `[doi:...](https://doi.org/...)` link, optionally followed by a replication-repo link (see below). Entries appear in reverse-chronological order within each subheading. Forthcoming entries use the pattern `*Journal*. Forthcoming.` with no DOI.
- **Verifying DOIs:** OUP pages (academic.oup.com) block direct fetches with HTTP 403. Resolve/verify book and chapter DOIs through the Crossref API instead: `https://api.crossref.org/works/<DOI>`.
- **Reading-list syllabi (e.g., `syllabi/SOCIOL204.qmd`):** Top-level `##` theme sections contain `###` group subheadings (e.g., "Concepts of Social Structure", "Cognition and Dual-Process Models"), each followed by a flat `-` bulleted list of APA-style citations. When adding or restructuring readings, group them thematically under new `###` subheadings and keep every citation verbatim.

## Replication repo links (2026-11-13)

- **CV format:** Papers with a public replication repo get `[{{< fa brands github >}} replication](https://github.com/olizardo/<repo>)` appended after the DOI link inside the entry paragraph. Use `{{< fa brands github >}}` (brands set) — plain `{{< fa github >}}` renders `fa-solid` and shows a tofu box for brand icons.
- **Finding the repo:** Match CV titles to `github.com/olizardo` repos by checking each repo's README (fetch `https://raw.githubusercontent.com/olizardo/<repo>/main/README.md`). Don't guess from repo names alone; verify HTTP 200 before linking.
- **Homepage nav:** The "Reproducible Research" dropdown in `_quarto.yml` should list all repos linked in the CV (keep in sync when adding CV replication links).
- **Pending (2026-11-16):** The CV links `olizardo/childress-lizardo-aesthetic-politics` for the Childress & Lizardo "An Impure Gaze" (*Socius*, forthcoming) entry, but the repo returned HTTP 404 at write time (not yet public). Recheck before touching that entry; create or rename the repo if the link should change.
- **Icon overlap bug (fixed 2026-11-13, do not regress):** FA's default `display: inline-block` on `<i>` collapses the glyph's advance width under `.cventry p { text-indent: -2.5rem }` in Chrome. `styles.css` forces `.cventry p i.fa-brands/-solid/-regular { display: inline }` — keep that rule.

## PDF/FontAwesome toolchain (2026-11-11)

- Quarto's fontawesome extension emits `fontawesome6` for LaTeX; system TeX Live 2022 lacks it, so FontAwesome 6 was installed manually from CTAN into the user tree (`~/texmf`, with `updmap-user --enable Map=fontawesome6.map`). If PDF rendering fails with `fontawesome6.sty not found`, re-install from `https://mirrors.ctan.org/fonts/fontawesome6.zip` (copy tex/, tfm/, type1/, map/, enc/, opentype/, doc/ into `~/texmf`, then `mktexlsr ~/texmf` + updmap-user).
- cv.qmd renders to both HTML and PDF (`_site/cv.html`, `_site/cv.pdf`); verify both after render, and check `pdffonts _site/cv.pdf` shows FontAwesome6 when icons are used.
- Verify rendered CV/pages visually with chrome-headless-shell screenshots (see user memory for the command); talls screenshots (`--window-size=1280,4600`) plus PIL crops work well for the CV's publication sections.

## Title capitalization convention (2026-11-16)

- Article, chapter, and talk titles in `cv.qmd` are sentence case: capitalize the first word and proper nouns only (CoDiNG, NetSense, NetHealth, Fitbit, Hawkes, Big Five, Erickson, etc.); after a colon, capitalize the first word of the subtitle. Verified and normalized across Publications, Invited Talks, and other quoted titles on 2026-11-16 — keep new entries in this style.
- Deliberately kept capitalized: stylized system/acronym names (SoCS, DHB, HNDS-R, "Physical Activity Trend eXtraction"), the idiom "The Ol' Ball and Chain", "Comment on <Name>" descriptors, and titles of books/reviewed works (which stay title case in italics).
- Organized-session titles in Professional Service follow sentence case too ("Recent advances in mathematical sociology"), except proper ASA section names in parentheses; session title "The making and remaking of Science..." deliberately keeps capital-S "Science" (discipline-stylization, flagged as ambiguous 2026-11-16).
