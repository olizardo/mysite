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

## DOI coverage in cv.qmd (2026-11-16)

- Every entry in the DOI-capable venues of the publication sections (Books, Peer-Reviewed Articles, Book Chapters, Comments & Short Pieces, Encyclopedia Entries) now carries a `[doi:...](https://doi.org/...)` link; verified via the Crossref API on 2026-11-16. Forthcoming entries intentionally have no DOI.
- DOIs added this session: Sage encyclopedia entries (Habitus 10.4135/9781452276052.n152; Cultural capital 10.4135/9781412994248.n144), "Foreign terrorism" (Routledge reissue of Ciment's Social Issues in America, 10.4324/9781315700724-157), "Jean Piaget: Sociology beyond holism and individualism" (10.5040/9798216016236.ch-033, ABC-CLIO/Bloomsbury AHFE-style chapter DOI), "Impact of network connectivity and agent commitment..." (10.54941/ahfe100197 — registered under AHFE International though printed by CRC Press), "Taste and the logic of practice in Distinction" (10.13060/00380288.2014.50.3.105, replacing the JSTOR link), and Weininger & Lizardo "Introduction" in Ritual, Emotion, Violence (10.4324/9780429464157-1).
- Verified DOI-less (do not keep hunting): Bergesen & Lizardo "Terrorism and world system theory" (2002 PISM volume, print-only), "Types of terrorism by world system location" (Humboldt Journal of Social Relations, JSTOR only), Sepulvado et al. trajectory-similarity chapter, "Embodied culture as procedure" (Warde & Southerton volume), the "Mapping the architectonic" review, and all newsletter/foreword/HBR-letter/AAAI-proceedings/Connections items (venues without DOIs).
- Crossref API rate-limits (HTTP 429) under parallel webfetch bursts; space queries a few seconds apart (e.g., curl + sleep in a loop) when checking many entries.

## Session updates (2026-11-20)

- New CV replication link: Lizardo & Skiles, "Musical taste and patterns of symbolic exclusion in the United States 1993–2012" (*Poetics* 53) → `olizardo/lizardo-skiles-age-field-musical-taste` (matched by README, not repo name; `lizardo-skiles-categorical-tolerance` is the *Sociological Science* paper, `skiles-lizardo-taste-influence-paper` is a different under-review manuscript). Reproducible Research dropdown in `_quarto.yml` updated to match.

## Publish workflow (2026-11-20)

- The site is deployed by `.github/workflows/publish.yml` ("Publish to GitHub Pages"): on every push to `main` it checks out the repo, installs Quarto + TinyTeX, renders the full site, and deploys to GitHub Pages. Concurrency group `pages` cancels superseded runs, so rapid pushes are safe.
- **No local full-site renders needed for deployment:** push to `main` and let CI build. Per-file `quarto render <file>.qmd` is still the right way to check a single document locally before pushing; a full local `quarto render` is optional and slow (CV renders both HTML and PDF).
