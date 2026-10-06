-- Wrap .cventry divs in the cventry LaTeX environment (PDF output only)
function Div(el)
  if el.classes:includes('cventry') and FORMAT:match('latex') then
    local blocks = el.content
    blocks:insert(1, pandoc.RawBlock('tex', '\\begin{cventry}'))
    blocks:insert(pandoc.RawBlock('tex', '\\end{cventry}'))
    return blocks
  end
end
