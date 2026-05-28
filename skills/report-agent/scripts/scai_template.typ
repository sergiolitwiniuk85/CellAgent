// SCAI Report Template — Single Cell AI Framework
// Usage: pandoc report.md -o report.typ --to typst --template scai_template.typ
//        typst compile report.typ

#let scai-title(project-name, date, session-id) = {
  align(center)[
    #block(text(weight: "bold", size: 2.5em)[#project-name])
    #v(0.5em)
    #text(size: 1.2em, fill: gray)[Single-Cell Analysis Report]
    #v(0.3em)
    #text(size: 0.9em, fill: gray)[#date \ SCAI Session: #session-id]
    #v(0.5em)
    #line(length: 60%, stroke: 1pt)
  ]
}

#let scai-section(title) = {
  heading(1, outlined: true, numbering: "1")[#title]
}

#let scai-subsection(title) = {
  heading(2, outlined: true, numbering: "1.1")[#title]
}

#let scai-table(headers, rows) = {
  table(
    columns: (auto, auto),
    inset: 8pt,
    stroke: 0.5pt,
    table.header(
      for h in headers [
        *#h*
      ]
    ),
    for row in rows {
      for cell in row [#cell]
    }
  )
}

#let scai-metric-table(data) = {
  table(
    columns: (auto, auto, auto),
    inset: 6pt,
    stroke: 0.3pt,
    table.header(
      [*Metric*], [*Value*], [*Notes*]
    ),
    for (metric, value, note) in data {
      [ #metric ], [ #value ], [ #note ]
    }
  )
}

#let scai-note(body) = {
  block(
    fill: luma(240),
    inset: 8pt,
    radius: 4pt,
    stroke: 0.5pt + blue,
  )[#body]
}

#let scai-code-block(code) = {
  block(
    fill: luma(245),
    inset: 8pt,
    radius: 4pt,
    stroke: 0.3pt + luma(200),
    width: 100%,
  )[#text(font: "Liberation Mono", size: 0.8em, code)]
}

#set page(
  margin: (left: 2.5cm, right: 2.5cm, top: 2cm, bottom: 2cm),
  numbering: "1",
  number-align: center,
)

#set heading(numbering: "1.1")
#set text(font: "Liberation Sans", size: 11pt)
#set par(justify: true, leading: 0.65em)
#show heading.where(level: 1): it => {
  v(1em)
  text(weight: "bold", size: 1.5em, fill: navy)[#it.body]
  v(0.3em)
  line(length: 100%, stroke: 1pt + navy)
  v(0.5em)
}
#show heading.where(level: 2): it => {
  v(0.8em)
  text(weight: "bold", size: 1.2em)[#it.body]
  v(0.3em)
}
