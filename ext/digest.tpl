\documentclass[]{article}
\usepackage[paperwidth=8.2in, paperheight=200in, voffset=0.1in, top=0.1in, left=0.1in, right=0.1in]{geometry}
\usepackage{fontspec}
\setmonofont{DejaVu Sans Mono}   % Set the monospaced font
\usepackage{newunicodechar}
\usepackage{amssymb,amsmath}
\usepackage{bm}
\usepackage{graphicx}
\usepackage{microtype}
\usepackage{hyperref}
{{#:tex_deps}}
{{{ :tex_deps }}}
{{/:tex_deps}}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1.2ex}
\usepackage{listings}
\lstset{
  columns=fullflexible,
  keepspaces=true,
  basicstyle=\ttfamily,
  breaklines=false
}

\hypersetup
       {   pdfauthor = { {{{:author}}} },
           pdftitle={ {{{:title}}} },
           colorlinks=TRUE,
           linkcolor=black,
           citecolor=blue,
           urlcolor=blue
       }

{{#:title}}
\title{ {{{ :title }}} }
{{/:title}}

{{#:author}}
\author{ {{{ :author }}} }
{{/:author}}

{{#:date}}
\date{ {{{ :date }}} }
{{/:date}}

{{ :highlight }}

\begin{document}

{{#:title}}\maketitle{{/:title}}

\vskip -0.5cm

{{{ :body }}}

\end{document}
