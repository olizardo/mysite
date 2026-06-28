# R script to build cv.qmd from LaTeX sources
library(utils)

# Source directory for CV files
src_dir <- "vita-src"

# Define the text cleaning function
clean_tex_to_md <- function(filepath, is_list = FALSE) {
  if (!file.exists(filepath)) {
    warning(paste("File does not exist:", filepath))
    return("")
  }
  
  lines <- readLines(filepath, warn = FALSE)
  # Remove lines that are purely LaTeX comments
  lines <- lines[!grepl("^\\s*%", lines)]
  text <- paste(lines, collapse = "\n")
  
  # Remove inline comments (?<!\\)%[^\n]*
  text <- gsub("(?<!\\\\)%[^\n]*", "", text, perl = TRUE)
  
  # 1. Normalize typos and line-splits in source LaTeX files before any parsing!
  text <- gsub("uclablue\\}\\s*\\n\\s*", "uclablue\\}", text, perl = TRUE)
  text <- gsub("\\bxtcolor\\{", "\\\\textcolor\\{", text, perl = TRUE)
  text <- gsub("\\begin{sloppypar}", "", text, fixed = TRUE)
  text <- gsub("\\end{sloppypar}", "", text, fixed = TRUE)
  
  # Normalize math superscripts inside textcolor arguments so they match footnotetext markers exactly!
  text <- gsub("{$^{\\dag}$}", "{\\dag}", text, fixed = TRUE)
  text <- gsub("{$^{\\dag \\dag}$}", "{\\dag \\dag}", text, fixed = TRUE)
  text <- gsub("{$^{\\dag \\dag \\dag}$}", "{\\dag \\dag \\dag}", text, fixed = TRUE)
  text <- gsub("{$^{\\dag \\dag \\dag \\dag}$}", "{\\dag \\dag \\dag \\dag}", text, fixed = TRUE)
  text <- gsub("{$^\\dag$}", "{\\dag}", text, fixed = TRUE)
  text <- gsub("{$^{\\ddag}$}", "{\\ddag}", text, fixed = TRUE)
  
  # 2. Remove layout/spacing/formatting commands completely at the very beginning
  # so they don't corrupt brace-matching for bold and italics!
  text <- gsub("(?m)^[ \\t]*\\\\setlength.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\itemsep.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\multicolsep.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\begin\\{multicols\\}.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\end\\{multicols\\}.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\begin\\{itemize\\}.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\end\\{itemize\\}.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\vspace\\{[^}]+\\}$", "", text, perl = TRUE)
  text <- gsub("\\\\vspace\\{[^}]+\\}", "", text, perl = TRUE) # inline vspace
  text <- gsub("(?m)^[ \\t]*\\\\footnotesize.*$", "", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]*\\\\small.*$", "", text, perl = TRUE)
  text <- gsub("\\\\footnotesize", "", text, fixed = TRUE)
  text <- gsub("\\\\small", "", text, fixed = TRUE)
  
  # 3. Specific footnotetext / mismatched markers replacements BEFORE general replacements!
  text <- gsub("\\footnotetext{\\textcolor{uclablue}{*}Co-Chair.}", "* *Co-Chair.*", text, fixed = TRUE)
  text <- gsub("\\footnotetext{\\textcolor{uclablue}{*}Honors Thesis.}", "* *Honors Thesis.*", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{\\ddag}\\footnotetext{\\textcolor{uclablue}{\\ddag}Theology.}", "^[Theology.]", text, fixed = TRUE)
  
  # 4. General footnote conversion using exact backreference matches
  pat <- "\\\\textcolor\\{([a-zA-Z]+)\\}\\{(.*?)\\}\\\\footnotetext\\{\\\\textcolor\\{\\1\\}\\{\\2\\}\\s*(.*?)\\}"
  text <- gsub(pat, "^[\\3]", text, perl = TRUE)
  
  # 5. Clean up superscript commas between multiple footnotes
  text <- gsub("$^{,}$", "", text, fixed = TRUE)
  
  # 6. Handle standard/residual standalone textcolor markers (not followed by footnotetext)
  text <- gsub("\\textcolor{uclablue}{*}", "*", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{**}", "**", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{\\ddag}", "‡", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{\\dag}", "<sup>†</sup>", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{\\dag \\dag}", "<sup>††</sup>", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{\\dag \\dag \\dag}", "<sup>†††</sup>", text, fixed = TRUE)
  text <- gsub("\\textcolor{uclablue}{\\dag \\dag \\dag \\dag}", "<sup>††††</sup>", text, fixed = TRUE)
  
  # 7. List items - must run BEFORE dash replacements!
  text <- gsub("\\item[--]", "-", text, fixed = TRUE)
  text <- gsub("\\item[- ]", "-", text, fixed = TRUE)
  text <- gsub("\\item[-]", "-", text, fixed = TRUE)
  text <- gsub("\\item", "-", text, fixed = TRUE)
  
  # 8. Replace newline with space (to avoid unrendered/escaped LaTeX command and overfull margins!)
  text <- gsub("\\newline", " ", text, fixed = TRUE)
  
  # 9. Clean up indentation of bullet points and text to prevent them from becoming monospace code blocks!
  # Use [ \t] instead of \s so we do NOT match or consume newline characters \n!
  # This completely preserves double newlines \n\n (empty lines) between paragraph blocks!
  text <- gsub("(?m)^[ \\t]+-\\s+", "  - ", text, perl = TRUE)
  text <- gsub("(?m)^[ \\t]+(?!-)", "", text, perl = TRUE)
  text <- gsub("(?m)^\\t+-\\s+", "  - ", text, perl = TRUE)
  text <- gsub("\t", "  ", text, fixed = TRUE)
  
  # 10. Simple fixed replacements for LaTeX commands/accents
  text <- gsub("\\&", "&", text, fixed = TRUE)
  text <- gsub("\\_", "_", text, fixed = TRUE)
  text <- gsub("\\l{}", "ł", text, fixed = TRUE)
  text <- gsub("\\l ", "ł ", text, fixed = TRUE)
  text <- gsub("\\l", "ł", text, fixed = TRUE)
  text <- gsub("\\'{a}", "á", text, fixed = TRUE)
  text <- gsub("\\'a", "á", text, fixed = TRUE)
  text <- gsub("\\'{o}", "ó", text, fixed = TRUE)
  text <- gsub("\\'o", "ó", text, fixed = TRUE)
  text <- gsub("\\\"{o}", "ö", text, fixed = TRUE)
  text <- gsub("\\\"o", "ö", text, fixed = TRUE)
  text <- gsub("\\v{s}", "š", text, fixed = TRUE)
  text <- gsub("\\v s", "š", text, fixed = TRUE)
  text <- gsub("\\`{a}", "à", text, fixed = TRUE)
  text <- gsub("\\`a", "à", text, fixed = TRUE)
  text <- gsub("\\\"{u}", "ü", text, fixed = TRUE)
  text <- gsub("\\\"u", "ü", text, fixed = TRUE)
  text <- gsub("\\~{n}", "ñ", text, fixed = TRUE)
  text <- gsub("\\~n", "ñ", text, fixed = TRUE)
  text <- gsub("\\c{c}", "ç", text, fixed = TRUE)
  text <- gsub("\\c c", "ç", text, fixed = TRUE)
  text <- gsub("\\v{c}", "ć", text, fixed = TRUE)
  text <- gsub("\\'{c}", "ć", text, fixed = TRUE)
  text <- gsub("\\'{n}", "ń", text, fixed = TRUE)
  text <- gsub("\\o{}", "ø", text, fixed = TRUE)
  
  # Typography
  text <- gsub("---", "—", text, fixed = TRUE)
  text <- gsub("--", "–", text, fixed = TRUE)
  text <- gsub("``", "\"", text, fixed = TRUE)
  text <- gsub("''", "\"", text, fixed = TRUE)
  text <- gsub("`", "'", text, fixed = TRUE)
  
  # LaTeX layout spacing - convert \ind to double newlines to split paragraphs!
  text <- gsub("\\noindent", "", text, fixed = TRUE)
  text <- gsub("\\medskip", "\n", text, fixed = TRUE)
  text <- gsub("\\bigskip", "\n\n", text, fixed = TRUE)
  text <- gsub("\\newpage", "", text, fixed = TRUE)
  text <- gsub("\\ind ", "\n\n", text, fixed = TRUE)
  text <- gsub("\\ind", "\n\n", text, fixed = TRUE)
  
  # Dynamic regex matching for bold, italic, href
  text <- gsub("\\\\textbf\\{([^}]+)\\}", "**\\1**", text, perl = TRUE)
  text <- gsub("\\\\emph\\{([^}]+)\\}", "*\\1*", text, perl = TRUE)
  text <- gsub("\\\\textit\\{([^}]+)\\}", "*\\1*", text, perl = TRUE)
  text <- gsub("\\\\nolinkurl\\{([^}]+)\\}", "\\1", text, perl = TRUE)
  
  # Handle the dual backslash or braces for bold/italic LaTeX structures
  text <- gsub("\\{\\\\bf\\s+([^}]+)\\}", "**\\1**", text, perl = TRUE)
  text <- gsub("\\{\\\\em\\s+([^}]+)\\}", "*\\1*", text, perl = TRUE)
  text <- gsub("\\{\\\\it\\s+([^}]+)\\}", "*\\1*", text, perl = TRUE)
  
  # Match \href{url}{text}
  text <- gsub("\\\\href\\{([^}]+)\\s*\\}\\{([^}]+)\\s*\\}", "[\\2](\\1)", text, perl = TRUE)
  
  # Remove remaining double backslashes \\ (LaTeX newlines) with empty string
  text <- gsub("\\\\", "", text, fixed = TRUE)
  
  # 11. Convert raw text list items ending in \\ or standalone plain text lines in list-blocks
  # into standard Markdown bullet list items! This runs on the clean Markdown text.
  if (is_list) {
    text <- gsub("(?m)^[ \\t]*([A-Za-z][^\\n]*?)$", "- \\1", text, perl = TRUE)
  }
  
  # Clean up multiple newlines (collapse 3 or more into exactly 2)
  text <- gsub("\\n{3,}", "\n\n", text)
  
  return(trimws(text))
}

# Start writing Quarto document
qmd_content <- c()

# Add Quarto YAML Front Matter
qmd_content <- c(qmd_content, "---")
qmd_content <- c(qmd_content, "title: \"Curriculum Vitae\"")
qmd_content <- c(qmd_content, "subtitle: \"Omar Lizardo\"")
qmd_content <- c(qmd_content, "format:")
qmd_content <- c(qmd_content, "  html:")
qmd_content <- c(qmd_content, "    toc: true")
qmd_content <- c(qmd_content, "    toc-depth: 3")
qmd_content <- c(qmd_content, "    css: styles.css")
qmd_content <- c(qmd_content, "  pdf:")
qmd_content <- c(qmd_content, "    pdf-engine: xelatex")
qmd_content <- c(qmd_content, "    documentclass: article")
qmd_content <- c(qmd_content, "    geometry: \"margin=1in\"")
qmd_content <- c(qmd_content, "    include-in-header:")
qmd_content <- c(qmd_content, "      - text: |")
qmd_content <- c(qmd_content, "          \\usepackage{parskip}")
qmd_content <- c(qmd_content, "          \\usepackage{fancyhdr}")
qmd_content <- c(qmd_content, "          \\usepackage{hanging}")
qmd_content <- c(qmd_content, "          \\usepackage{xcolor}")
qmd_content <- c(qmd_content, "          \\usepackage{xurl}")
qmd_content <- c(qmd_content, "          \\usepackage{microtype}")
qmd_content <- c(qmd_content, "          \\usepackage{hyperref}")
qmd_content <- c(qmd_content, "          \\hypersetup{breaklinks=true, colorlinks=true, urlcolor=ucladark}")
qmd_content <- c(qmd_content, "          \\definecolor{uclablue}{HTML}{005587}")
qmd_content <- c(qmd_content, "          \\definecolor{ucladark}{HTML}{003B5C}")
qmd_content <- c(qmd_content, "          \\definecolor{uclagold}{HTML}{FFB81C}")
qmd_content <- c(qmd_content, "          \\definecolor{uclamagenta}{HTML}{FF00A5}")
qmd_content <- c(qmd_content, "          \\pagestyle{fancy}")
qmd_content <- c(qmd_content, "          \\fancyhf{}")
qmd_content <- c(qmd_content, "          \\renewcommand{\\headrulewidth}{0pt}")
qmd_content <- c(qmd_content, "          \\rhead{\\footnotesize Last updated: \\today}")
qmd_content <- c(qmd_content, "          \\lfoot{\\footnotesize Omar Lizardo — Curriculum Vitae}")
qmd_content <- c(qmd_content, "          \\cfoot{\\footnotesize \\thepage}")
qmd_content <- c(qmd_content, "          \\newenvironment{cventry}{\\begin{hangparas}{2.5em}{1}}{\\end{hangparas}}")
qmd_content <- c(qmd_content, "---")
qmd_content <- c(qmd_content, "")

# Add elegant HTML download button and hide it in PDF
qmd_content <- c(qmd_content, "::: {.content-visible when-format=\"html\"}")
qmd_content <- c(qmd_content, "[{{< fa file-pdf >}} Download PDF Version](vita.pdf){.btn .btn-outline-primary .btn-sm role=\"button\"}")
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")

# Add CV contact header
qmd_content <- c(qmd_content, "::: {.text-center .content-visible when-format=\"pdf\"}")
qmd_content <- c(qmd_content, "# Omar Lizardo")
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")

qmd_content <- c(qmd_content, "::: {.small .text-muted .text-center .mb-4}")
qmd_content <- c(qmd_content, "Department of Sociology, University of California, Los Angeles  ")
qmd_content <- c(qmd_content, "264 Haines Hall, 375 Portola Plaza, Los Angeles, CA, 90095  ")
qmd_content <- c(qmd_content, "[olizardo@soc.ucla.edu](mailto:olizardo@soc.ucla.edu) | [olizardo.bol.ucla.edu](http://olizardo.bol.ucla.edu) | [GitHub](https://github.com/olizardo) | [Google Scholar](https://scholar.google.com/citations?user=Lt5nMNkAAAAJ) | [ORCID: 0000-0002-5405-3007](https://orcid.org/0000-0002-5405-3007)")
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")

# 1. Education
qmd_content <- c(qmd_content, "## Education")
qmd_content <- c(qmd_content, "::: {.cventry}")
qmd_content <- c(qmd_content, "*University of Arizona (Tucson, AZ)*")
qmd_content <- c(qmd_content, "**PhD, Sociology**, 2006.")
qmd_content <- c(qmd_content, "**MA, Sociology**, 2002.")
qmd_content <- c(qmd_content, "")
qmd_content <- c(qmd_content, "*Brooklyn College, City University of New York (Brooklyn, NY)*")
qmd_content <- c(qmd_content, "**BS, Psychology**, 1997.")
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")

# Helper to add a hanging-indent section
add_hanging_section <- function(title, filename) {
  content <- clean_tex_to_md(file.path(src_dir, filename))
  if (content != "") {
    res <- c(paste("##", title), "::: {.cventry}", content, ":::", "")
    return(res)
  }
  return(c())
}

# Helper to add a standard (list-based) section
add_standard_section <- function(title, filename) {
  content <- clean_tex_to_md(file.path(src_dir, filename))
  if (content != "") {
    res <- c(paste("##", title), content, "")
    return(res)
  }
  return(c())
}

# 2. Academic Positions
qmd_content <- c(qmd_content, add_hanging_section("Academic Positions", "01-positions.tex"))

# 3. Professional Memberships
qmd_content <- c(qmd_content, add_standard_section("Professional Memberships", "03-memberships.tex"))

# 4. Distinctions, Honors, and Awards
qmd_content <- c(qmd_content, add_standard_section("Distinctions, Honors, and Awards", "02-awards.tex"))

# 5. Books
qmd_content <- c(qmd_content, add_hanging_section("Books", "04-books.tex"))

# 6. Peer-Reviewed Articles
qmd_content <- c(qmd_content, "## Peer-Reviewed Articles")
# We list files in order
peer_reviewed_files <- list(
  list(title = "In Press / Recent", file = "05-current.tex"),
  list(title = "2024–2025", file = "05-2024-2025.tex"),
  list(title = "2022–2023", file = "05-2022-2023.tex"),
  list(title = "2021", file = "05-2021.tex"),
  list(title = "2019–2020", file = "05-2019-2020.tex"),
  list(title = "2018", file = "05-2018.tex"),
  list(title = "2017", file = "05-2017.tex"),
  list(title = "2016", file = "05-2016.tex"),
  list(title = "2015", file = "05-2015.tex"),
  list(title = "2013–2014", file = "05-2013-2014.tex"),
  list(title = "2010–2012", file = "05-2010-2012.tex"),
  list(title = "2009", file = "05-2009.tex"),
  list(title = "2008", file = "05-2008.tex"),
  list(title = "2007", file = "05-2007.tex"),
  list(title = "2006", file = "05-2006.tex"),
  list(title = "2003–2005", file = "05-2003-2005.tex")
)

for (item in peer_reviewed_files) {
  content <- clean_tex_to_md(file.path(src_dir, item$file))
  if (content != "") {
    qmd_content <- c(qmd_content, paste("###", item$title), "::: {.cventry}", content, ":::", "")
  }
}

# 7. Book Chapters
qmd_content <- c(qmd_content, "## Book Chapters")
chapter_files <- list(
  list(title = "2021–2030", file = "06-chapters-2021-2030.tex"),
  list(title = "2016–2020", file = "06-chapters-2016-2020.tex"),
  list(title = "2002–2015", file = "06-chapters-2002-2015.tex")
)
for (item in chapter_files) {
  content <- clean_tex_to_md(file.path(src_dir, item$file))
  if (content != "") {
    qmd_content <- c(qmd_content, paste("###", item$title), "::: {.cventry}", content, ":::", "")
  }
}

# 8. Comments & Short Pieces
qmd_content <- c(qmd_content, "## Comments & Short Pieces")
comment_files <- list(
  list(title = "2016–Present", file = "07-comments-2016-present.tex"),
  list(title = "2011–2014", file = "07-comments-2011-2014.tex"),
  list(title = "2008–2010", file = "07-comments-2008-2010.tex")
)
for (item in comment_files) {
  content <- clean_tex_to_md(file.path(src_dir, item$file))
  if (content != "") {
    qmd_content <- c(qmd_content, paste("###", item$title), "::: {.cventry}", content, ":::", "")
  }
}

# 9. Encyclopedia Entries
qmd_content <- c(qmd_content, add_hanging_section("Encyclopedia Entries", "08-encyclopedia.tex"))

# 10. Book Reviews
qmd_content <- c(qmd_content, "## Book Reviews")
review_files <- list(
  list(title = "2020–Present", file = "09a-reviews-2020-present.tex"),
  list(title = "2011–2019", file = "09b-reviews-2011-2019.tex"),
  list(title = "2003–2010", file = "09c-reviews-2003-2010.tex")
)
for (item in review_files) {
  content <- clean_tex_to_md(file.path(src_dir, item$file))
  if (content != "") {
    qmd_content <- c(qmd_content, paste("###", item$title), "::: {.cventry}", content, ":::", "")
  }
}

# 11. arXiv / SocArXiv e-prints
qmd_content <- c(qmd_content, add_hanging_section("arXiv e-prints", "10a-arxiv.tex"))
qmd_content <- c(qmd_content, add_hanging_section("SocArXiv e-prints", "10b-socarxiv.tex"))

# 12. Invited Lectures and Addresses
qmd_content <- c(qmd_content, "## Invited Lectures and Addresses")
lecture_files <- list(
  list(title = "2021–Present", file = "11-lec-2021-present.tex"),
  list(title = "2018–2020", file = "11-lec-2018-2020.tex"),
  list(title = "2015–2017", file = "11-lec-2015-2017.tex"),
  list(title = "2010–2014", file = "11-lec-2010-2014.tex"),
  list(title = "2004–2009", file = "11-lec-2004-2009.tex")
)
for (item in lecture_files) {
  content <- clean_tex_to_md(file.path(src_dir, item$file))
  if (content != "") {
    qmd_content <- c(qmd_content, paste("###", item$title), "::: {.cventry}", content, ":::", "")
  }
}

# 13. Grants and Sponsored Programs
qmd_content <- c(qmd_content, add_hanging_section("Grants and Sponsored Programs", "12-grants.tex"))

# 14. Conference Presentations
qmd_content <- c(qmd_content, "## Conference Presentations")
conf_files <- list(
  list(title = "2020–Present", file = "13-conf-2020-present.tex"),
  list(title = "2018–2019", file = "13-conf-2018-2019.tex"),
  list(title = "2016–2017", file = "13-conf-2016-2017.tex"),
  list(title = "2010–2015", file = "13-conf-2010-2015.tex"),
  list(title = "2006–2009", file = "13-conf-2006-2009.tex"),
  list(title = "2002–2005", file = "13-conf-2002-2005.tex")
)
for (item in conf_files) {
  content <- clean_tex_to_md(file.path(src_dir, item$file))
  if (content != "") {
    qmd_content <- c(qmd_content, paste("###", item$title), "::: {.cventry}", content, ":::", "")
  }
}

# 15. Editorial Work
qmd_content <- c(qmd_content, "## Editorial Work")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "14-editor.tex")), "")
qmd_content <- c(qmd_content, "### Editorial Boards")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "15-edboard.tex")), "")
qmd_content <- c(qmd_content, "### Program Committees")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "16-prog-comm.tex")), "")

# 16. Professional Service
qmd_content <- c(qmd_content, "## Professional Service")
qmd_content <- c(qmd_content, "### ASA Committees")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "17-asa-committee.tex")), "")
qmd_content <- c(qmd_content, "### Conference Organizer")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "18-organizer.tex")), "")
qmd_content <- c(qmd_content, "### Conference Session Discussant")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "19-discussant.tex")), "")
qmd_content <- c(qmd_content, "### Conference Session Presider")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "20-presider.tex")), "")

# 17. Student Advising
qmd_content <- c(qmd_content, "## Student Advising")
qmd_content <- c(qmd_content, "### PhD Dissertation Committee")
qmd_content <- c(qmd_content, "::: {.cventry}")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "21a-phd.tex"), is_list = TRUE))
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")
qmd_content <- c(qmd_content, "### MA Thesis Committee")
qmd_content <- c(qmd_content, "::: {.cventry}")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "21b-ma.tex"), is_list = TRUE))
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")
qmd_content <- c(qmd_content, "### Undergraduate Honors Thesis Advising")
qmd_content <- c(qmd_content, "::: {.cventry}")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "21c-ug.tex"), is_list = TRUE))
qmd_content <- c(qmd_content, ":::")
qmd_content <- c(qmd_content, "")

# 18. Courses Taught
qmd_content <- c(qmd_content, "## Courses Taught")
qmd_content <- c(qmd_content, clean_tex_to_md(file.path(src_dir, "22-courses.tex"), is_list = TRUE), "")

# Write out the file
writeLines(qmd_content, "cv.qmd")
cat("SUCCESS: cv.qmd has been built!\n")
