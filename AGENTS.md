# Agent Instructions

This is the repository for my academic website, built using Quarto.
When interacting with this codebase, please adhere to the following rules:

1. **Quarto:** Remember that this site relies on Quarto. Any changes to navigation should be made in `_quarto.yml`. 
2. **Re-rendering:** Run `quarto render` or `quarto render <file/folder>` to build and test changes to the website structure or content.
3. **Data/Content Migrations:** Keep consistent formatting for migrated content (e.g. keeping YAML frontmatter clean and minimal, using the standard `.qmd` extensions).
4. **CV Generation:** The CV is a single, self-contained Quarto document. Any updates, additions, or modifications to the CV should be made directly in `cv.qmd` using standard Markdown formatting.
5. **Themes/CSS:** For aesthetic or visual adjustments, utilize the custom stylesheets (`custom.scss`, `custom-dark.scss`, and `styles.css`).

## Content conventions (2026-11-03)

- **cv.qmd:** Publications are organized as `##` sections (Books, Peer-Reviewed Articles, Book Chapters, Comments & Short Pieces, Encyclopedia Entries, etc.) with `###` period subheadings; each entry is one paragraph ending in a `[doi:...](https://doi.org/...)` link. Entries appear in reverse-chronological order within each subheading.
- **Verifying DOIs:** OUP pages (academic.oup.com) block direct fetches with HTTP 403. Resolve/verify book and chapter DOIs through the Crossref API instead: `https://api.crossref.org/works/<DOI>`.
- **Reading-list syllabi (e.g., `syllabi/SOCIOL204.qmd`):** Top-level `##` theme sections contain `###` group subheadings (e.g., "Concepts of Social Structure", "Cognition and Dual-Process Models"), each followed by a flat `-` bulleted list of APA-style citations. When adding or restructuring readings, group them thematically under new `###` subheadings and keep every citation verbatim.
