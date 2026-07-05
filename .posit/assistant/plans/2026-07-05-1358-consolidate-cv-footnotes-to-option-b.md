# Plan - Consolidate CV Footnotes to Option B (Traditional Symbols)

This plan details how to consolidate the repetitive footnotes in `cv.qmd` ("Graduate Seminar", "Co-chair", "Computer Science and Engineering", and "Honors Thesis") into clean, reusable superscript symbols, removing 24 duplicate footnotes from the bottom of the document.

## 1. Symbol Mapping
We will use the following traditional footnote symbols:
* `^†^` (Dagger) — **Graduate Seminar** (replaces `[^gs]` and `^[Graduate Seminar.]`)
* `^*^` (Asterisk) — **Co-chair** (replaces `[^co]`)
* `^‡^` (Double dagger) — **Computer Science and Engineering** (replaces `[^cs]`)
* `^§^` (Section sign) — **Honors Thesis** (replaces `[^ho]`)

## 2. Step-by-Step Changes in `cv.qmd`

### Section: Student Advising (Lines 1332–1477)
* **Co-chair (`[^co]` -> `^*^`):**
  - Line 1343: `- Stephanie Zhang[^co] (2028)` -> `- Stephanie Zhang^*^ (2028)`
  - Line 1344: `- Nida Sanglimsuwan[^co] (2028)` -> `- Nida Sanglimsuwan^*^ (2028)`
  - Line 1349: `- Justin Farrell[^co] (2013)` -> `- Justin Farrell^*^ (2013)`
  - Line 1354: `- Michael Wood[^co] (2019)` -> `- Michael Wood^*^ (2019)`
  - Line 1355: `- Brandon Sepulvado[^co] (2019)` -> `- Brandon Sepulvado^*^ (2019)`
  - Line 1356: `- Dustin Stoltz[^co] (2020)` -> `- Dustin Stoltz^*^ (2020)`
  - Line 1358: Remove the old footnote definition `[^co]: Co-chair.`

* **Computer Science and Engineering (`[^cs]` -> `^‡^`):**
  - Line 1390: `- Yang Yang[^cs] (2015)` -> `- Yang Yang^‡^ (2015)`
  - Line 1394: `- Yuxiao Dong[^cs]  (2017)` -> `- Yuxiao Dong^‡^  (2017)`
  - Line 1396: `- Aastha Nigam[^cs]  (2018)` -> `- Aastha Nigam^‡^  (2018)`
  - Line 1401: Remove the old footnote definition `[^cs]: Computer Science and Engineering.`

* **Honors Thesis (`[^ho]` -> `^§^`):**
  - Line 1466: `- Teresa Bone[^ho] (2008)` -> `- Teresa Bone^§^ (2008)`
  - Line 1467: `- Melissa Truitt[^ho] (2011)` -> `- Melissa Truitt^§^ (2011)`
  - Line 1468: `- Yo Tam Yoon[^ho] (2012)` -> `- Yo Tam Yoon^§^ (2012)`
  - Line 1469: `- Olevia Boykin[^ho] (2014)` -> `- Olevia Boykin^§^ (2014)`
  - Line 1470: `- Karyn Vilbig[^ho] (2014)` -> `- Karyn Vilbig^§^ (2014)`
  - Line 1471: `- Molly Feeney[^ho] (2015)` -> `- Molly Feeney^§^ (2015)`
  - Line 1474: Remove the old footnote definition `[^ho]: Honors Thesis.`

* **Advising Section Legends:**
  At the end of the **Student Advising** section (around line 1474), we will add a clean grouped legend:
  ```markdown
  ^*^ Co-chair.  
  ^‡^ Computer Science and Engineering.  
  ^§^ Honors Thesis.
  ```

---

### Section: Courses Taught (Lines 1477–1499)
* **Graduate Seminar (`[^gs]` & `^[Graduate Seminar.]` -> `^†^`):**
  - Line 1480: `- Cultural Sociology^[Graduate Seminar.]` -> `- Cultural Sociology^†^`
  - Line 1481: `- Topics in Sociological Theorizing[^gs]` -> `- Topics in Sociological Theorizing^†^`
  - Line 1483: `- Social Network Methods (A & B)[^gs]` -> `- Social Network Methods (A & B)^†^`
  - Line 1484: `- Theory and Research in Sociology[^gs]` -> `- Theory and Research in Sociology^†^`
  - Line 1488: `- Classical Sociological Theory[^gs]` -> `- Classical Sociological Theory^†^`
  - Line 1490: `- Contemporary Theory[^gs]` -> `- Contemporary Theory^†^`
  - Line 1491: `- Culture and Cognition[^gs]` -> `- Culture and Cognition^†^`
  - Line 1494: `- From Publishable to Published[^gs]` -> `- From Publishable to Published^†^`
  - Line 1496: `- Social Networks[^gs]` -> `- Social Networks^†^`
  - Line 1498: Remove old footnote definition `[^gs]: Graduate Seminar.`

* **Teaching Section Legend:**
  At the bottom of the **Courses Taught** section (around line 1498), we will place:
  ```markdown
  ^†^ Graduate Seminar.
  ```

## 3. Build & Verification
We will run `quarto render cv.qmd` to verify:
1. The site compiles with zero errors or warnings.
2. The HTML/PDF layout renders beautifully.
3. The 24 duplicate footnotes are successfully consolidated and no longer pollute the bottom of the CV.
