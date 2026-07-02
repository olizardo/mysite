# Agent Instructions

This is the repository for my academic website, built using Quarto.
When interacting with this codebase, please adhere to the following rules:

1. **Quarto:** Remember that this site relies on Quarto. Any changes to navigation should be made in `_quarto.yml`. 
2. **Re-rendering:** Run `quarto render` or `quarto render <file/folder>` to build and test changes to the website structure or content.
3. **Data/Content Migrations:** Keep consistent formatting for migrated content (e.g. keeping YAML frontmatter clean and minimal, using the standard `.qmd` extensions).
4. **CV Generation:** Changes to my CV should be made in the LaTeX files within `vita-src/`, and then built by running `Rscript build_cv.R`.
5. **Themes/CSS:** For aesthetic or visual adjustments, utilize the custom stylesheets (`custom.scss`, `custom-dark.scss`, and `styles.css`).
