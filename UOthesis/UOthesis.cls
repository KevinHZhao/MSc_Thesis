%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the file UOthesis.cls
%
% Inspired by  MPLETTER.CLS  from Guide to Latex, 4th edition
% by H. Kopka and P. W. Daly
%
% Copyright 2008 by Benoit Dionne
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}[1996/06/01]
\ProvidesClass{UOthesis}[2023/01/01 v0.003 UO thesis style]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Packages required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{ifthen}
\RequirePackage{amsmath}
\RequirePackage{amsthm}
\RequirePackage{amssymb}
\RequirePackage{amscd}        % For simple commutative diagrams
\RequirePackage{graphicx}
\RequirePackage{rotating}     % To rotate figures, tables, ...
\RequirePackage{color}
\RequirePackage{pdfpages}     % To include research papers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Possible options: all normal options for report style plus
% french, Nomencl, Makeindex, NoLofF and NoLofT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newboolean{@french}
\setboolean{@french}{false}
\DeclareOption{french}{\setboolean{@french}{true}}

\newboolean{@Nomencl}      %  To include the list of symbols
\setboolean{@Nomencl}{false}
\DeclareOption{Nomencl}{\setboolean{@Nomencl}{true}}

\newboolean{@Makeindex}      %  To include an index at the end of the thesis
\setboolean{@Makeindex}{false}
\DeclareOption{Makeindex}{\setboolean{@Makeindex}{true}}

\newboolean{@NoLofF}       %  To exclude the list of figures
\setboolean{@NoLofF}{false}
\DeclareOption{NoLofF}{\setboolean{@NoLofF}{true}}

\newboolean{@NoLofT}       %  To exclude the list of tables
\setboolean{@NoLofT}{false}
\DeclareOption{NoLofT}{\setboolean{@NoLofT}{true}}

% The option should be used to get proper page headers.
% It is passed to the report class by default.
\DeclareOption*{\PassOptionsToClass{\CurrentOption,twoside}{report}}
\ProcessOptions
\LoadClass{report}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Packages required for the PDF/A documents.
% The result is still not a fully compliant PDF/A document.
% The newest version of pdftex seems to produce by default what the
% following packages are supposed to do.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
\RequirePackage{filecontents}
\RequirePackage[utf8x]{inputenc}
\RequirePackage{colorprofiles}
\RequirePackage[T1]{fontenc}
\RequirePackage[a-3b]{pdfx}
%\RequirePackage[a-2b,mathxmp]{pdfx}[2018/12/22]

% hyperref is included by pdfx but pdfx doesn't complain if hyperref
% is not installed on the system.  It simply doesn't do anything.
% This line will produce the complaint if hyperref is not installed.
\RequirePackage{hyperref}
\hypersetup{pdfencoding=unicode}
% \hypersetup{pdfstartview=}
% KZ Adding cleveref for easier referencing figs/tables, using \cref
% instead of \ref
\usepackage[nameinlink,capitalize]{cleveref}
%% A space is still included when multiple sections are referenced,
%% but I think removing the space for multiple section refs would look
%% worse.
\crefname{chapter}{Chapter}{Chapters}
\Crefname{chapter}{Chapter}{Chapters}
\crefname{section}{\S}{\S\S}
\Crefname{section}{Section}{Sections}
\crefname{equation}{Equation}{Equations}
\Crefname{equation}{Equation}{Equations}
\crefname{figure}{Figure}{Figures}
\Crefname{figure}{Figure}{Figures}
\Crefname{appsec}{Appendix}{appendices}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure of the theorems, propositions, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \theoremstyle{plain}
\theoremstyle{definition}
\ifthenelse{\boolean{@french}}{
  \newtheorem{theo}{Th\'eor\`eme}[section]
  \newtheorem{lem}[theo]{Lemme}
  \newtheorem{defn}[theo]{D\'efinition}
  \newtheorem{cor}[theo]{Corollaire}
  \newtheorem{prop}[theo]{Proposition}
}{
  \newtheorem{theo}{Theorem}[section]
  \newtheorem{lem}[theo]{Lemma}
  \newtheorem{defn}[theo]{Definition}
  \newtheorem{cor}[theo]{Corollary}
  \newtheorem{prop}[theo]{Proposition}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Counter for the equation number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setcounter{secnumdepth}{3}
\numberwithin{equation}{section}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To allow page breaks in multi-line equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\allowdisplaybreaks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package to create a list of symbols if one is to be created.
% The command
%    \nomenclature[prefix]{symbol}{description⟩}
% must be given after the definition of each symbol.
% prefix is optional and use as a prefix before sorting.  It does not
% appear in the list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{@Nomencl}}{
  \RequirePackage[refpage,intoc]{nomencl}
  \makenomenclature
}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package to create the index.  The command
%    \index{main_entry[!sub_entry][!subsub_entry]}
% must be inserted after the main_entry[!sub_entry][!subsub_entry]
% that you want to index.   The text between [...] is optional.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{@Makeindex}}{
  \RequirePackage{makeidx}
  \makeindex
}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table of contents, list of figures, list of tables, and list of
% symbols.
% If your thesis has section numbers larger than 9, many tables or
% many figures, a lot of pages, ... then some of the items on a line
% of the table of contents, the list of tables, the list of figures or
% the list of symbols may overlap.  The following commands should
% solve these problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newlength{\TofCvalueA}
\setlength{\TofCvalueA}{2em}
\renewcommand\l@section{\@dottedtocline{1}{\TofCvalueA}{3em}}  % space
                                % for the section number of the table
                                % of contents
\newcommand{\addSpaceContents}[1]{\addtolength{\TofCvalueA}{#1}}

\newlength{\TofCvalueB}
\newlength{\TofCvalueC}
\setlength{\TofCvalueB}{3em}
\setlength{\TofCvalueC}{4em}
\renewcommand\@pnumwidth{\TofCvalueB} % space for the page number of
                                % the table of contents
\renewcommand\@tocrmarg{\TofCvalueC}  % space for the right margin of
                                % the table of contents
\newcommand{\addSpacePageNumber}[1]{
  \addtolength{\TofCvalueB}{#1}
  \addtolength{\TofCvalueC}{#1}
}
\newcommand{\addSpaceRightMargin}[1]{
  \addtolength{\TofCvalueC}{#1}
  \ifthenelse{\TofCvalueB<\TofCvalueC}{
    \setlength{\TofCvalueB}{\TofCvalueC}
  }{}
}

\ifthenelse{\boolean{@Nomencl}}{
  \setlength{\nomitemsep}{-\parsep}
  \newlength{\NomenCvalueA}
  \setlength{\NomenCvalueA}{1.5cm}
  \renewcommand{\nomlabelwidth}{\NomenCvalueA}  % space for the label of the
                                % symbol in a list of symbols 
  \newcommand{\addSpaceLabelSymbols}[1]{\addtolength{\NomenCvalueA}{#1}}

  \newlength{\NomenCvalueB}
  \setlength{\NomenCvalueB}{1.5cm}
  \def\Dotfill{\xleaders\hbox to1em{\hfill.\hfill}\hfill}  % space
                                % between the dots of a dotted line
  \renewcommand{\pagedeclaration}[1]{\ \Dotfill \makebox[\NomenCvalueB][r]{#1}}
                                % format for the page number of the list
  \newcommand{\addSpacePageSymbols}[1]{\addtolength{\NomenCvalueB}{#1}}
}{}

\ifthenelse{\boolean{@NoLofF}}{}{
  \newlength{\LofFvalue}
  \setlength{\LofFvalue}{1.5em}
  \newcommand{\addSpaceFigure}[1]{
    \addtolength{\LofFvalue}{#1}
  }
  \renewcommand\l@figure[2]{\@dottedtocline{1}{1.5em}{\LofFvalue}{#1}{#2}}
                                      % space for the figure number
}

\ifthenelse{\boolean{@NoLofT}}{}{
  \newlength{\LofTvalue}
  \setlength{\LofTvalue}{1.5em}
  \newcommand{\addSpaceTable}[1]{
    \addtolength{\LofTvalue}{#1}
  }
  \renewcommand\l@table[2]{\@dottedtocline{1}{1.5em}{\LofTvalue}{#1}{#2}}
                                      % space for the table number
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Il est possible de compiler un texte qui possede des accents, ...
% sans avoir a recourir aux symboles \'e, \"u, etc.
% En plus de permettre l'utilisation d'un editeur de texte normal,
% vous pouvez aussi utiliser votre "spell checker" prefere.
%
% Table of contents, list of figures, list of tables, and list of
% symbols, ...  in French and English.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{@french}}{
  \RequirePackage[french]{babel}
  \RequirePackage{ucs}
  % \RequirePackage[T1]{fontenc}       % Already included
  % \RequirePackage[utf8x]{inputenc}   % Already included
  % \RequirePackage[latin1]{inputenc}  % Requis selon votre
                                       % environment de travail

  \ifthenelse{\boolean{@Nomencl}}{
    \renewcommand{\nomname}{Liste des symboles}
  }{}
  \renewcommand{\contentsname}{Table des mati\`eres}
  \renewcommand{\listfigurename}{Liste des figures}
  \renewcommand{\listtablename}{Liste des tableaux}
  \newcommand{\dedicationsname}{D\'edicaces}
  \newcommand{\acknowledgementsname}{Remerciements}
  \newcommand{\prefacename}{Pr\'eface}
}{
  \ifthenelse{\boolean{@Nomencl}}{
    \renewcommand{\nomname}{List of Symbols}
  }{}
  \renewcommand{\contentsname}{Contents}
  \renewcommand{\listfigurename}{List of Figures}
  \renewcommand{\listtablename}{List of Tables}
  \newcommand{\dedicationsname}{Dedications}
  \newcommand{\acknowledgementsname}{Acknowledgements}
  \newcommand{\prefacename}{Preface}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Page style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\raggedbottom
\RequirePackage[letterpaper,lmargin=1in,height=9in,width=6.5in,top=1in]{geometry}
\setlength\parskip{1ex plus 0.3ex minus 0ex}
\setlength\parindent{2em}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% font for the fancy header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\hdf}{\sffamily\bfseries}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Header
%
% 1) the page number on the left followed by the title of the
%    chapter for the even numbered pages (left pages).
% 2) the page number on the right preceded by the title of the
%    section for the odd numbered pages (right pages).
% For a one-sided document, only 1 is applied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\setlength{\headheight}{1.5\baselineskip}
\setlength{\headsep}{1.5\baselineskip}
\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{\thechapter.\ #1}}{}}
\renewcommand{\sectionmark}[1]{\markright{\MakeUppercase{\thesection.\ #1}}}
\renewcommand{\headrulewidth}{0.5pt}
\renewcommand{\plainheadrulewidth}{0pt}
\fancyhf{}
\fancyhead[LE,RO]{\hdf \thepage}
\fancyhead[LO]{\hdf \rightmark}
\fancyhead[RE]{\hdf \leftmark}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New environment for chapter without number
%
% \nonumchapter{chapter_name}
% to produce the header of a chapter without a number;
% the chapter is still listed in the table of contents.
%
% \notachapter{chapter_name}
% to produce the header of a chapter without a number;
% the chapter is not listed in the table of contents.
%
% Note that the numbering of equations, ... in a chapter with no
% chapter number will be based on the previous chapter and section
% numbers.  This is certainly not what we should have.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
\newcommand{\nonumchapter}[1]{
  \chapter*{#1 \@mkboth{\MakeUppercase{#1}}{\MakeUppercase{#1}}}
  \addcontentsline{toc}{chapter}{#1}
}

\newcommand{\notachapter}[1]{
  \chapter*{#1 \@mkboth{\MakeUppercase{#1}}{\MakeUppercase{#1}}}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important variables supplied by the user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\UO}{University of Ottawa}
\newcommand{\setUO}[1]{\renewcommand{\UO}{#1}}

\newcommand{\UOprogrfac}{Faculty of Science}
\newcommand{\setUOprogrfac}[1]{\renewcommand{\UOfac}{#1}}

\newcommand{\UOdept}{Department of Mathematics and Statistics}
\newcommand{\setUOdept}[1]{\renewcommand{\UOdept}{#1}}

\newcommand{\UOname}{}       % Full name of the candidate
\newcommand{\setUOname}[1]{\renewcommand{\UOname}{#1}}

\newcommand{\UOfdegree}{}    % Full title of the degree sought
\newcommand{\setUOfdegree}[1]{ \renewcommand{\UOfdegree}{#1} }

\newcommand{\UOdegree}{}     % Abbreviation of the title of the degree sought
\newcommand{\setUOdegree}[1]{ \renewcommand{\UOdegree}{#1} }

\newcommand{\UOtitlefootnote}{}   % title for the institute footnote
\newcommand{\setUOtitlefootnote}[1]{\renewcommand{\UOtitlefootnote}{#1}}

\newcommand{\UOcpryear}{}  % Year that should appear with the copyright notice
\newcommand{\setUOcpryear}[1]{\renewcommand{\UOcpryear}{#1}}

\newcommand{\UOtitle}{}    % The title of the thesis
\newcommand{\setUOtitle}[1]{\renewcommand{\UOtitle}{#1}}

% If French option has been chosen
\ifthenelse{\boolean{@french}}{
  \renewcommand{\UO}{Universit\'e d'Ottawa}
  \renewcommand{\UOgradfac}{Facult\'e des sciences}
  \renewcommand{\UOprogrfac}{Facult\'e des sciences}
  \renewcommand{\UOdept}{D\'epartement de math\'ematiques et de statistique}
}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short cut to define all the variables associated to the degree
% sought.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newboolean{undergrad}
\setboolean{undergrad}{false}
\newcommand{\UOgoal}{}
\newcommand{\setUOgoal}[1]{\renewcommand{\UOgoal}{#1}}
\newcommand{\phd}{
  \ifthenelse{\boolean{@french}}{
    \setUOfdegree{Doctorat en philosophie Math\'ematiques et statistique}
    \setUOtitlefootnote{doctorat}
  }{
    \setUOfdegree{Doctorate in Philosophy Mathematics and Statistics}
    \setUOtitlefootnote{Ph.D.}
  }
  \setUOdegree{Ph.D.}
  \ifthenelse{\boolean{@french}}{
    \setUOgoal{Th\`ese soumisse en vue de l'obtention du\\ \UOfdegree}
  }{
    \setUOgoal{Thesis submitted in partial fulfillment of the requirements for the degree of\\ \UOfdegree}
  }
}
\newcommand{\msc}{
  \ifthenelse{\boolean{@french}}{
    \setUOfdegree{Ma\^{\i}trise \`es sciences Math\'ematiques et statistique}
    \setUOtitlefootnote{ma\^{\i}trise}
  }{
    \setUOfdegree{Master of Science Mathematics and Statistics}
    \setUOtitlefootnote{M.Sc.}
  }
  \setUOdegree{M.Sc.}
  \ifthenelse{\boolean{@french}}{
    \setUOgoal{Th\`ese soumisse en vue de l'obtention de la\\ \UOfdegree}
  }{
    \setUOgoal{Thesis submitted in partial fulfillment of the requirements for the degree of\\ \UOfdegree}
  }
}
\newcommand{\project}{
  \ifthenelse{\boolean{@french}}{
    \setUOfdegree{Ma\^{\i}trise \`es sciences Math\'ematiques et statistique}
    \setUOtitlefootnote{ma\^{\i}trise}
  }{
    \setUOfdegree{Master of Science Mathematics and Statistics}
    \setUOtitlefootnote{M.Sc.}
  }
  \setUOdegree{M.Sc.}
  \ifthenelse{\boolean{@french}}{
    \setUOgoal{Projet soumis en vue de l'obtention de la\\ \UOfdegree}
  }{
    \setUOgoal{Project submitted in partial fulfillment of the requirements for the degree of\\ \UOfdegree}
  }
}
\newcommand{\uproject}{
  \ifthenelse{\boolean{@french}}{
    \setUOfdegree{Baccalaur\'eat}
    \setUOtitlefootnote{}
  }{
    \setUOfdegree{Bachelor}
    \setUOtitlefootnote{}
  }
  \setUOdegree{B.Sc.}
  \ifthenelse{\boolean{@french}}{
    \setUOgoal{Projet soumis en vue de l'obtention du \UOfdegree}
  }{
    \setUOgoal{Project submitted in partial fulfillment of the requirements of the \UOfdegree}
  }
  \setboolean{undergrad}{true}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The title page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\UOtitlepage}{
  \thispagestyle{empty}
  \vspace*{0.3in}
  \begin{sloppypar}
    \begin{center}
      \Large\upshape{\UOtitle}
    \end{center}
  \end{sloppypar}
  \vfill
  \begin{center}
    \large\rm \UOname
  \end{center}
  \vfill
  \begin{center}
    \renewcommand{\thefootnote}{\fnsymbol{footnote}}
    \ifthenelse{\boolean{@french}}{
      \UOgoal\ifthenelse{\boolean{undergrad}}{}{
        \footnote[1]{le programme de \UOtitlefootnote\ est un
          programme conjoint avec l'Universit\'e Carleton,
          administr\'e par l'Institut en math\'ematiques et en
          statistique d'Ottawa-Carleton}} \
    }{
      \UOgoal\ifthenelse{\boolean{undergrad}}{}
      {\footnote[1]{The \UOtitlefootnote\ program is a joint program
          with Carleton University, administered by the
          Ottawa-Carleton Institute of Mathematics and Statistics}}
    }
    \renewcommand{\thefootnote}{\arabic{footnote}}
  \end{center}
  \vfill
  \begin{center}
    \UOdept\\
    \UOprogrfac\\
    \UO
  \end{center}
  \vfill
  \begin{center}
    \copyright\ \UOname, Ottawa, Canada, \UOcpryear
  \end{center}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The bibliography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hand made bibliography
\newenvironment{BasicBibliography}{
  \begin{thebibliography}{999}
    \ifthenelse{\boolean{@french}}{
      \addcontentsline{toc}{chapter}{Bibliographie}
    }{
      \addcontentsline{toc}{chapter}{Bibliography}
    }
  }{\end{thebibliography}}

% Bibliography made with BibTeX
% Include only the references in the bib file which are cited in the
% text.
\newcommand{\nostarBibTex}[1]{
  \bibliography{#1}
  \ifthenelse{\boolean{@french}}{
    \addcontentsline{toc}{chapter}{Bibliographie}
  }{
    \addcontentsline{toc}{chapter}{Bibliography}
  }
}

% Bibliography made with BibTeX
% Include all references in the bib file even if they are not cited in
% the text.
\newcommand{\starBibTex}[1]{
  \nocite{*}
  \bibliography{#1}
  \ifthenelse{\boolean{@french}}{
    \addcontentsline{toc}{chapter}{Bibliographie}
  }{
    \addcontentsline{toc}{chapter}{Bibliography}
  }
}

\makeatletter
\newcommand{\bibTexBibliography}{\@ifstar\starBibTex\nostarBibTex}
\makeatother

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The proof environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The filled box corresponding to the empty box suggested by the AMS
% \renewcommand{\qedsymbol}{$\blacksquare$}
% I prefer this box
\renewcommand{\qedsymbol}{\rule{0.4em}{0.7em}}

% I prefer bold characters.
\renewcommand*{\proofname}{{\normalfont\bfseries Proof}}

\makeatletter
\renewenvironment{proof}[1][\proofname]{\par\pushQED{\qed}%
\normalfont \topsep6\p@\@plus6\p@\relax
\trivlist
\item\relax {\normalfont\bfseries #1\@addpunct{.}}%
\hspace\labelsep\ignorespaces
}{\popQED\endtrivlist\@endpefalse}
\makeatother

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commands to be run after \begin{document}
% The command \AtBeginDocument{} is not really after
% the \begin{document}.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
\AtBeginDocument{
  \pagenumbering{roman}
  \UOtitlepage
  \clearpage

  \IfFileExists{abstract}{\notachapter{Abstract}
    \input{abstract}
    \clearpage}{}
  \IfFileExists{resume}{\notachapter{R\'esum\'e}
    \input{resume}
    \clearpage}{}
  \IfFileExists{dedications}{\notachapter{\dedicationsname}
    \input{resume}
    \clearpage}{}
  \IfFileExists{acknowledgement}{\notachapter{\acknowledgementsname}
    \input{resume}
    \clearpage}{}

  \tableofcontents
  \clearpage

  \ifthenelse{\boolean{@NoLofF}}{}{
    \listoffigures
    \addcontentsline{toc}{chapter}{\listfigurename}
  }
  \ifthenelse{\boolean{@NoLofT}}{}{
    \listoftables
    \addcontentsline{toc}{chapter}{\listtablename}
  }
  \ifthenelse{\boolean{@Nomencl}}{
    \markboth{\MakeUppercase{\nomname}}{\MakeUppercase{\nomname}}
    \printnomenclature
  }{}
  \clearpage

  \IfFileExists{preface}{\notachapter{\prefacename}
    \input{preface}
    \clearpage}{}

  \pagenumbering{arabic}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commands to be run just before \begin{document}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\AtEndDocument{
  \ifthenelse{\boolean{@Makeindex}}{
    \addcontentsline{toc}{chapter}{\indexname}
    \printindex}{}
}
