%% latex thesis template
%% Copyright (C) 2015      Simon Woertz <firstname@lastname.at>
%% Copyright (C) 2009-2011 Christoph Sommer <christoph.sommer@uibk.ac.at>
%%                         Falko Dressler <falko.dressler@uibk.ac.at>


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{thesis-style}


% -----------------------------------------------------------------------------
% Class options
% -----------------------------------------------------------------------------
\RequirePackage{ifthen}

% [germanthesis] - Thesis is written in German
\newboolean{germanthesis}
\setboolean{germanthesis}{false}
\DeclareOption{germanthesis}{\setboolean{germanthesis}{true}}

% [plainunnumbered] - Don't print numbers on plain pages
\newboolean{plainunnumbered}
\setboolean{plainunnumbered}{false}
\DeclareOption{plainunnumbered}{\setboolean{plainunnumbered}{true}}

% [earlydraft] - Settings for quick draft printouts
\newboolean{earlydraft}
\setboolean{earlydraft}{false}
\DeclareOption{earlydraft}{\setboolean{earlydraft}{true}}

% [watermark] - Print current time/date at bottom of each page
\newboolean{watermark}
\setboolean{watermark}{false}
\DeclareOption{watermark}{\setboolean{watermark}{true}}

% [phdthesis] - switch to PhD thesis style
\newboolean{phdthesis}
\setboolean{phdthesis}{false}
\DeclareOption{phdthesis}{\setboolean{phdthesis}{true}}

% [twoside] - double sided
\newboolean{oneside}
\setboolean{oneside}{true}
\DeclareOption{twoside}{\setboolean{oneside}{false}}

% [cutmargins] - text body fills complete page
\newboolean{cutmargins}
\setboolean{cutmargins}{false}
\DeclareOption{cutmargins}{\setboolean{cutmargins}{true}}

\ProcessOptions \relax
% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Basics
% -----------------------------------------------------------------------------
% This class is based on report
\ifthenelse{\boolean{oneside}}{
\LoadClass[a4paper,11pt,onecolumn,oneside,openright]{report}
}{
\LoadClass[a4paper,11pt,onecolumn,twoside,openright]{report}
}
\ifthenelse{\boolean{cutmargins}}{
\RequirePackage[paperwidth=365pt,paperheight=618pt,textwidth=345pt,textheight=598pt,bindingoffset=0cm,marginparwidth=3cm,hmarginratio=1:1,vmarginratio=1:1]{geometry}
}{
\RequirePackage[a4paper,textwidth=345pt,textheight=598pt,bindingoffset=0cm,marginparwidth=3cm,hmarginratio=1:1,vmarginratio=1:1]{geometry}
}

% 150 percent line spacing
\RequirePackage{setspace}
\onehalfspacing

% Allow more (and larger) floats on text pages
\renewcommand{\topfraction}{0.9}  % use up to ..% of space on top
\renewcommand{\bottomfraction}{0.8}  % use up to ..% of space on bottom
\setcounter{topnumber}{4} % place up to .. on top
\setcounter{bottomnumber}{1} % place up to .. on bottom
\setcounter{totalnumber}{4} % place up to .. total
\renewcommand{\textfraction}{0.07} % allow down to ..% of text
\renewcommand{\floatpagefraction}{0.7} % fill at least ..% of float pages (must be less than \topfraction)

% Avoid widows and orphans (single line on bottom or top of page, respectively) at almost any cost
\clubpenalty = 10000
\widowpenalty = 10000
%\raggedbottom

% Default fonts
\newcommand{\mypageheadfont}{\normalfont}
\newcommand{\myheadingfont}{\normalfont}
%\RequirePackage[osf]{libertine}  % Linux Libertine
%\renewcommand{\mypageheadfont}{\normalfont\libertine}
%\renewcommand{\myheadingfont}{\normalfont\libertine}
\RequirePackage{lmodern}  % TT: Latin Modern
\RequirePackage{helvet}  % SF: Helvetica 
\RequirePackage[bitstream-charter,sfscaled=false]{mathdesign}  % RM: Bitstream Charter
% TODO: nicer math font?

% More encoding and typesetting fixes and tweaks
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{textcomp}
\RequirePackage{fixltx2e}
\ifthenelse{\boolean{earlydraft}}{
	\newcommand{\mydisableprotrusion}{}
}{
	\RequirePackage[babel,protrusion=true,expansion=true,tracking=false,kerning=true,spacing=true]{microtype}
	\newcommand{\mydisableprotrusion}{\microtypesetup{protrusion=false}}
}
\RequirePackage{hyphenat}

% Graphics insertion and filetypes' order of precedence
\RequirePackage[pdftex]{graphicx}
\pdfminorversion=6
\DeclareGraphicsRule{.out.png}{png}{.out.png}{}
\ifthenelse{\boolean{earlydraft}}{
	\DeclareGraphicsExtensions{.thumb.jpg,.out.png,.pdf,.png,.jpg}
	\DeclareGraphicsRule{.thumb.jpg}{jpg}{.thumb.jpg}{}
}{
	\DeclareGraphicsExtensions{.out.png,.pdf,.png,.jpg}
}

% Typesetting URLs
\RequirePackage[lowtilde,hyphens]{url}
\urlstyle{tt}

% For hyperref<->algorithm interaction
\RequirePackage{float}

% For function listing
\RequirePackage{titletoc}

% Cross-references are clickable
\RequirePackage[pdftex,bookmarks=true,bookmarksnumbered=true,colorlinks=false,pdfdisplaydoctitle=true,breaklinks=true]{hyperref}

% Captions
\RequirePackage[font=small, labelfont=bf, labelsep=endash, margin=1cm]{caption}

% Subfigures
%\RequirePackage[margin=0cm]{subfig}

% Citations
%\RequirePackage{cite}
\RequirePackage[backend=bibtex, style=numeric-comp]{biblatex}

% Quotes
% German
\RequirePackage[autostyle]{csquotes}

% The "SI" and "num" commands for typesetting units
\RequirePackage[range-phrase={\dots\,},per-mode=symbol,detect-all,load-configurations=binary,forbid-literal-units]{siunitx}
\ifthenelse{\boolean{germanthesis}}{
	\sisetup{locale=DE}
}{
}
\DeclareSIUnit\byte{Byte}
\DeclareSIUnit\decibelm{dBm}
\DeclareSIUnit\mph{mph}

% Commands for beautiful tables
\RequirePackage{booktabs}

% Equation typesetting fixes and tweaks + proof/theorem environment
\RequirePackage{amsmath}
\RequirePackage{amsthm}
\newtheorem{thm}{Theorem}[section]

% Allow list iteration over comma separated parameters
\RequirePackage{pgf}
\RequirePackage{pgffor}

% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Internationalization support
% -----------------------------------------------------------------------------

% Multi-language support using Babel
\let\myoldshow\show\renewcommand{\show}{}  % temporarily disable \show command, as babeltools 2010/10/14 seems to use it
\RequirePackage[shorthands=off]{babeltools}  % don't use any shorthands (i.e. "a will not try to produce an umlaut)
\let\show\myoldshow\renewcommand{\myoldshow}{}  % re-enable \show command
\ifthenelse{\boolean{germanthesis}}{
	\RequirePackage[american,ngerman]{babel}
}{
	\RequirePackage[ngerman,american]{babel}
}

% % Don't indent paragraphs in german theses. Instead, leave blank space between paragraphs.
% \ifthenelse{\boolean{germanthesis}}{
% 	\RequirePackage{parskip}
% }{
% }

% Time and date functions
\RequirePackage[useregional]{datetime2}

% Clever References
\ifthenelse{\boolean{germanthesis}}{
	\RequirePackage[ngerman]{cleveref}[2009/17/04]
	\def\cref{\Cref}
	\def\crefrange{\Crefrange}
	\crefname{lstlisting}{Listing}{Listings}
	\Crefname{lstlisting}{Listing}{Listings}
	\crefname{listing}{Listing}{Listings}
	\Crefname{listing}{Listing}{Listings}
	\crefname{algorithm}{Algorithmus}{Algorithmen}
	\Crefname{algorithm}{Algorithmus}{Algorithmen}
}{
	\RequirePackage{cleveref}[2009/17/04]
	\def\cref{\Cref}
	\def\crefrange{\Crefrange}
	\crefname{lstlisting}{listing}{listings}
	\Crefname{lstlisting}{Listing}{Listings}
	\crefname{listing}{listing}{listings}
	\Crefname{listing}{Listing}{Listings}
	\crefname{algorithm}{algorithm}{algorithms}
	\Crefname{algorithm}{Algorithm}{Algorithms}
}

% Some more language-dependent commands
\ifthenelse{\boolean{germanthesis}}{
	%\newcommand{\thedate}{\selectlanguage{ngerman}\today\selectlanguage{american}}
	\newcommand{\thedate}{\today}
	\newcommand{\glossarytitlename}{Abkürzungsverzeichnis}
}{
	\newcommand{\thedate}{\today}
	\newcommand{\glossarytitlename}{List of Acronyms}
}

% Language-independent timestamp
\newcommand{\timestamp}{\renewcommand{\dateseparator}{-}\yyyymmdddate\today\quad\settimeformat{hhmmsstime}\currenttime}

% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% All the simple stuff goes here
% -----------------------------------------------------------------------------

% The Acronym package
\RequirePackage[printonlyused]{acronym}

% Make automatically inserted "blank pages" really blank
\RequirePackage{emptypage}

% Fix spacing around "..."
\RequirePackage{ellipsis}

% TeX ist kein Zeichenprogramm
\RequirePackage{tikz}[2010/10/13]
\usetikzlibrary{}

% Ensure that no floats are typeset before they are referenced
\RequirePackage{flafter}

% Provide the \FloatBarrier command that flushes all pending floats
\RequirePackage{placeins}

% Nag about deprecated LaTeX packages and commands
\RequirePackage[l2tabu, orthodox]{nag}

% Define command that will return current filename
\RequirePackage{currfile}

% Do not reset footnote counter on new chapter
\RequirePackage[]{chngcntr}
\counterwithout{footnote}{chapter}

% "Draft" watermark
\newcommand{\thewatermark}{\ifthenelse{\boolean{watermark}}{%
	\begin{tikzpicture}[remember picture, overlay]
		\node[anchor=south, outer sep=3.75cm] at (current page.south) {
			\color{gray}\currfilename\quad\timestamp
		};
	\end{tikzpicture}%
}{}}

% The "algorithm" environment
\RequirePackage{algorithmic}[2009/08/24]
%\RequirePackage[plain,chapter]{algorithm}
\RequirePackage[chapter]{algorithm}
\floatplacement{algorithm}{bp}

% The "syntax" environment
\RequirePackage[nounderscore,rounded]{syntax}

\renewcommand\fs@plain{%
\def\@fs@cfont{\rmfamily}%
\let\@fs@capt\floatc@plain%
\def\@fs@pre{\hrule\relax\kern2pt}%
\def\@fs@post{}%
\def\@fs@mid{\kern2pt\hrule\relax\vspace\abovecaptionskip\relax}%
\let\@fs@iftopcapt\iffalse%
}
\floatstyle{plain}
\restylefloat{algorithm}
\ifthenelse{\boolean{germanthesis}}{
	\floatname{algorithm}{Algorithmus}
}{
}

% Chapter/section headings and page styles
%\RequirePackage[pagestyles,nobottomtitles*]{titlesec}
\RequirePackage[pagestyles]{titlesec}
\ifthenelse{\boolean{plainunnumbered}}{
\renewpagestyle{plain}{%
	\sethead{}{}{}
	\setfoot{}{}{\thewatermark}
}
}{
\renewpagestyle{plain}{%
	\sethead{}{}{}
	\setfoot{}{\mypageheadfont\thepage}{\thewatermark}
}
}

% Pagestyle "fancy"
\newpagestyle{fancy}{%
\sethead*{
	\bfseries
	\toptitlemarks
	\mypageheadfont\bfseries
	\ifthechapter
		{
			\ifthesection
			{\thesection~\sectiontitle}
			{\thechapter~\chaptertitle}
		}
		{\chaptertitle}
}
{}{\mypageheadfont\bfseries\thepage}
	\setfoot{}{}{\thewatermark}
	\headrule
}
\pagestyle{fancy}

% rename (sub)paragraph to (sub)subsubsubsection
\let\myparagraph\paragraph
\renewcommand{\paragraph}{\FIXME{in order to avoid confusion with hyperref, \\paragraph has been renamed to \\subsubsubsection}}
\let\subsubsubsection\myparagraph
\let\mysubparagraph\subparagraph
\renewcommand{\subparagraph}{\FIXME{in order to avoid confusion with hyperref, \\subparagraph has been renamed to \\subsubsubsubsection}}
\let\subsubsubsubsection\mysubparagraph

% heading styles
\titleformat{\chapter}[display]{\setlength\parskip{0pt}\myheadingfont\huge\bfseries\titlerule\vspace{0.25cm}}{\LARGE\chaptertitlename~\thechapter}{1em}{}[\vspace{0.25cm}\titlerule]
\titleformat{\section}{\myheadingfont\Large\bfseries}{\thesection}{1em}{}
\titleformat{\subsection}{\myheadingfont\large\bfseries}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\myheadingfont\normalsize\bfseries}{\thesubsubsection}{1em}{}
\titleformat{\paragraph}{\myheadingfont\normalsize\bfseries}{\theparagraph}{1em}{}
\titleformat{\subparagraph}{\myheadingfont\normalsize}{\thesubparagraph}{1em}{\labelitemi~}
\titlespacing*{\chapter}{0pt}{50pt}{40pt minus \parskip}
\titlespacing*{\section}{0pt}{3.5ex plus 1ex minus .2ex - \parskip}{2.3ex plus .2ex - \parskip}
\titlespacing*{\subsection}{0pt}{3.25ex plus 1ex minus .2ex - \parskip}{1.5ex plus .2ex - \parskip}
\titlespacing*{\subsubsection}{0pt}{3.25ex plus 1ex minus .2ex - \parskip}{1.5ex plus .2ex - \parskip} 
\titlespacing*{\paragraph}{0pt}{3.25ex plus 1ex minus .2ex - \parskip}{1.5ex plus .2ex - \parskip} 
\titlespacing*{\subparagraph}{0pt}{3.25ex plus 1ex minus .2ex - \parskip}{1.5ex plus .2ex - \parskip} 
\ifthenelse{\boolean{earlydraft}}{
	\setcounter{tocdepth}{3}
	\setcounter{secnumdepth}{3}
}{
	\setcounter{tocdepth}{3}
	\setcounter{secnumdepth}{3}
}
% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Commands for setting author, title, ... 
% -----------------------------------------------------------------------------
\renewcommand{\author}[1]{\gdef\theauthor{#1}}
\renewcommand{\title}[1]{\gdef\thetitle{#1}}
\newcommand{\germantitle}[1]{\gdef\thegermantitle{#1}}
\newcommand{\thesistype}[1]{\gdef\thethesistype{#1}}
\newcommand{\thesistypegerman}[1]{\gdef\thethesistypegerman{#1}}
\newcommand{\thesiscite}[1]{\gdef\thethesiscite{#1}}
\newcommand{\birthday}[1]{\gdef\thebirthday{#1}}
\newcommand{\birthplace}[1]{\gdef\thebirthplace{#1}}
\newcommand{\thesisstart}[1]{\gdef\thethesisstart{#1}}
\newcommand{\thesisend}[1]{\gdef\thethesisend{#1}}
\newcommand{\advisors}[1]{\gdef\theadvisors{#1}}
\newcommand{\dean}[1]{\gdef\thedean{#1}}
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Listings
% -----------------------------------------------------------------------------
\RequirePackage[savemem]{listings}
\ifthenelse{\boolean{germanthesis}}{
	\renewcommand{\lstlistlistingname}{Listingverzeichnis}
}{
	\renewcommand{\lstlistlistingname}{List of Listings}
}
\lstset{
	language=C,
	frame=lines,
	framesep=5pt,
	captionpos=b,
	abovecaptionskip=1em,
	numbers=left,
	xleftmargin=15pt,
	framexleftmargin=15pt,
	numberstyle=\tiny,
	numbersep=5pt,
	stepnumber=1,
	fontadjust,
	tabsize=2,
	showtabs=false,
	showspaces=false,
	showstringspaces=false,
	breaklines=true,
	breakatwhitespace=true,
	keywordstyle=\color{black},
	identifierstyle=\bfseries\color{black},
	commentstyle=\color{blue},
	prebreak=\raisebox{0ex}[0ex][0ex]{\color{gray}\ensuremath{\searrow}},
	stringstyle=\bfseries\color{black},
	basicstyle=\ttfamily\bfseries\color{black}\footnotesize,
}
\lstdefinestyle{txt}{
	float,
	floatplacement=bp,
	language=C
}
\lstdefinestyle{xml}{
	float,
	floatplacement=bp,
	language=XML
}
% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% TOCs
% -----------------------------------------------------------------------------

% Include list-of-Xs in table of contents
\ifthenelse{\boolean{earlydraft}}{
	\RequirePackage[]{tocbibind}
	\renewcommand{\lstlistoflistings}{\begingroup\tocfile{\lstlistlistingname}{lol}\endgroup}
	\newcommand{\listofabbreviations}{\addcontentsline{toc}{chapter}{\glossarytitlename}}
	\newcommand{\mymaketodotocentry}{\todototoc}
}{
	\RequirePackage[nottoc,notlot,notlof]{tocbibind}
	\newcommand{\listofabbreviations}{}
	\newcommand{\mymaketodotocentry}{}
}

% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Visible TODO and FIXME markers
% -----------------------------------------------------------------------------
\ifthenelse{\boolean{germanthesis}}{
	\RequirePackage[colorinlistoftodos, textwidth=3.75cm, shadow, ngerman]{todonotes}
	\@todonotes@SetTodoListName{Offene Punkte}
}{
	\RequirePackage[colorinlistoftodos, textwidth=3.75cm, shadow]{todonotes}
}
\newcounter{TODOCOUNT}
\newcommand{\TODO}[1]{\vspace{0.5em}\todo[inline, color=orange]{#1}\stepcounter{TODOCOUNT}}
\newcommand{\FIXME}[1]{\todo[size=\small, color=red]{#1}\stepcounter{TODOCOUNT}}
\AtEndDocument{
	\ifnum\value{TODOCOUNT}>0
		%\cleardoublepage
		\mymaketodotocentry\listoftodos
	\fi
}
% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Title page
% -----------------------------------------------------------------------------
\newcommand{\coverpagetextsmall}[1]{%
	\@coverpage@smalltextsize%
	#1%
	\vspace{0.5cm}\\%
}

\newcommand{\coverpagetextlarge}[1]{%
	\vspace{0.25cm}%
	\@coverpage@largetextsize%
	#1%
	\@coverpage@smalltextsize%
	\vspace{0.5cm}\\%
}
\renewcommand{\maketitle}{%

% Do some stuff we had to postpone until after e.g. \thetitle was set
\hypersetup{%
  pdftitle = {\thetitle},
  pdfsubject = {\thethesistype},
  pdfauthor = {\theauthor},
  pdfborder=0 0 0
}

\pagenumbering{alph}

% Output cover page
\ifthenelse{\boolean{cutmargins}}{
	% no cover page if cutmargins is set
}
{%

\tikzoption{text @coverpage@ badly raggedleft}[]{\def\tikz@text@action{\raggedleft\relax}}

\begingroup
\pagestyle{empty}
\setcounter{page}\@ne

\def\@coverpage@titlesize{\@setfontsize{\@coverpage@titlesize}{20}{24pt}}
\def\@coverpage@subtitlesize{\@setfontsize{\@coverpage@subtitlesize}{14}{17pt}}
\def\@coverpage@citeassize{\@setfontsize{\@coverpage@citeassize}{8}{9pt}}
\def\@coverpage@addressize{\@setfontsize{\@coverpage@addressize}{8}{9pt}}
\def\@coverpage@smalltextsize{\@setfontsize{\@coverpage@smalltextsize}{14}{17pt}}
\def\@coverpage@largetextsize{\@setfontsize{\@coverpage@largetextsize}{20}{24pt}}

\renewcommand{\baselinestretch}{1.00}
\definecolor{uibk_darkgray}{HTML}{878787}
\definecolor{uibk_gray}{HTML}{DADADA}
\definecolor{uibk_green}{HTML}{ACC23C}
\begin{tikzpicture}[remember picture,overlay] \node[xshift=0cm,yshift=0cm] at (current page.south west) { \begin{tikzpicture}[remember picture, overlay]

	%%% logos
	\ifthenelse{\boolean{earlydraft}}{
		% no watermark in [earlydraft] mode
	}{
		\pgftext[bottom,right,at={\pgfpoint{25.5cm}{-5.5cm}}] {%
			\includegraphics[height=20cm]{logos/watermark_uibk}
		}
	}

	%%% Kaestchenleisten
	%\foreach \x in {0,1,...,17,20} {
	\foreach \x in {1.5,2.5,...,15.5,18.5} {
		\foreach \y in {23.5,24.5} {%
			\draw[fill=uibk_gray, draw=none] (0.05cm + \x * 1cm,\y * 1cm) rectangle +(0.9cm,0.9cm);
		}
	}
	\draw[fill=uibk_green, draw=none] (0.05cm + 2.5 * 1cm, 23.5 * 1cm) rectangle +(0.9cm,0.9cm);
	\pgftext[bottom,left,at={\pgfpoint{0.05cm + 16.5 * 1cm}{23.5 * 1cm}}] {%
		\includegraphics[width=1.9cm]{logos/logo_uibk_unofficial}%
	}
	\draw [uibk_gray, ultra thick] (1.5,4) -- (19.5,4) ;
	%\foreach \x in {1.5,2.0,...,19.0} {
	%	\foreach \y in {5} {%
	%	%\foreach \y in {4.5,19.6} {%
	%		\draw[fill=uibk_gray, draw=none] (0.05cm + \x * 1cm,\y * 1cm) rectangle +(0.4cm,0.4cm);
	%	}
	%}
	%\pgftext[top,left,at={\pgfpoint{2cm}{4cm}}] {%
	%	\includegraphics[height=2cm]{additional logo}%
	%}

	%%% single-line text
	\pgftext[top,left,at={\pgfpoint{2.5cm}{27cm}}] {%
		\raggedright%
		\sffamily%
		\bfseries%
		\@coverpage@subtitlesize%
		\renewcommand{\baselinestretch}{1.00}%
		\textcolor{uibk_darkgray}{universit\"{a}t innsbruck}%
	}
	\pgftext[bottom,left,at={\pgfpoint{2.5cm}{25.75cm}}] {%
		\raggedright%
		\sffamily%
		\bfseries%
		\@coverpage@subtitlesize%
		\renewcommand{\baselinestretch}{1.00}%
		\textcolor{uibk_green}{Fakult\"{a}t f\"{u}r Mathematik, Informatik und Physik}%
		%\textcolor{uibk_green}{Faculty of Mathematics, Computer Science and Physics}%
	}

	%%% multi-line text
	\pgftext[top,left,at={\pgfpoint{2.5cm}{20cm}}] {%
		\parbox{15cm}{%
			\raggedright%
			\sffamily%
			\@coverpage@largetextsize%
			\renewcommand{\baselinestretch}{1.25}%
			\coverpagetextsmall{\nohyphens{\theauthor}}%
			\coverpagetextlarge{\nohyphens{\thetitle}}%
			\coverpagetextsmall{\thethesistypegerman}%
			\coverpagetextsmall{\thedate}%
		}%
	};
%	\draw[xshift=2.5cm,yshift=5cm] node [above right,text width=15cm, text badly ragged] {%
%		\sffamily%
%		\@coverpage@citeassize%
%		\renewcommand{\baselinestretch}{1.00}%
%%		\begin{otherlanguage*}{ngerman}
%%			Bitte zitieren als:\\%
%%			\nohyphens{\theauthor, ``\thetitle,'' \thesistypegerman, Universit\"{a}t Innsbruck, Institut f\"{u}r Informatik, \monthname~\the\year.}%
%%		\end{otherlanguage*}
%		\begin{otherlanguage*}{american}%
%			Please cite as:\\%
%			\nohyphens{\theauthor, ``\thetitle,'' \thethesiscite, University of Innsbruck, Institute of Computer Science, \monthname~\the\year.}%
%		\end{otherlanguage*}
%	};
	\draw[xshift=19.08cm,yshift=3.5cm] node [below left,text width=12cm, text @coverpage@ badly raggedleft] {%
		\sffamily%
		\@coverpage@addressize%
		\renewcommand{\baselinestretch}{1.00}%
		%University of Innsbruck\\%
		%Institute of Computer Science\\%
		Universit\"{a}t Innsbruck\\%
		Institut f\"{u}r Informatik\\%
		%Computer and Communication Systems\\%
		\vspace{0.5em}%
		%Technikerstr. 21a $\cdot$ 6020 Innsbruck $\cdot$ Austria\\%
		Technikerstr. 21a $\cdot$ 6020 Innsbruck $\cdot$ \"{O}sterreich\\%
		\textbf{\textcolor{uibk_green}{\href{http://informatik.uibk.ac.at/}{http://informatik.uibk.ac.at/}}}\\%
	};
\end{tikzpicture} }; \end{tikzpicture}

\cleardoublepage
\setcounter{page}\@ne
\endgroup
}

% Title pages use single line spacing, no parindent and no parskips
\begin{singlespace}
\setlength\parindent{0pt} 
\setlength\parskip{0pt} 

%
% Now output the title page
%

\pagestyle{empty}
\pagenumbering{roman}

\ifthenelse{\boolean{phdthesis}}
{
\thispagestyle{empty}

%TODO: PhD thesis title page
}
% else: Master Thesis title page
{
\thispagestyle{empty}

\vspace*{0.8cm}

\begin{center}

{\LARGE \textbf{\thetitle}}\\[2\baselineskip]

{\large
\thethesistypegerman\\[2\baselineskip]

vorgelegt von\\[2\baselineskip]

\textbf{\theauthor}\\[1\baselineskip]

geb. am \thebirthday\\
in \thebirthplace\\[2\baselineskip]

angefertigt am\\[2\baselineskip]

\textbf{%Lehrstuhl für Technische Informatik\\[0.2\baselineskip]
%Computer and Communication Systems\\[1.2\baselineskip]
Institut für Informatik\\[0.2\baselineskip]
Leopold-Franzens-Universität Innsbruck\\[0.2\baselineskip]
}
}~\\[1\baselineskip]

\begin{tabular}{>{\raggedleft}p{5cm}p{6cm}}
Betreuer:	& %
  \foreach \n [count=\ni] in \theadvisors {%
    \ifnum\ni=1%
    \else%
      \newline
    \fi%
    \textbf{\n}%
  }%
\\
~ & ~\\
Abgabe der Arbeit: & \textbf{\selectlanguage{ngerman}\today\selectlanguage{american}}\\
\end{tabular}

\end{center}

\clearpage
}
\end{singlespace}

%
% The delaration comes here
%

\ifthenelse{\boolean{phdthesis}}{}
{
\thispagestyle{empty}

\begin{otherlanguage*}{ngerman}
{\Large\noindent%
Erklärung \newline}

\noindent%
Ich versichere, dass ich die Arbeit ohne fremde Hilfe und ohne Benutzung anderer als der angegebenen Quellen angefertigt habe und dass die Arbeit in gleicher oder ähnlicher Form noch keiner anderen Prüfungsbehörde vorgelegen hat und von dieser als Teil einer Prüfungsleistung angenommen wurde.\\
Alle Ausführungen, die wörtlich oder sinngemäß übernommen wurden, sind als solche gekennzeichnet.
\end{otherlanguage*}


\vspace{2cm}


\begin{otherlanguage*}{american}
{\Large\noindent%
Declaration \newline}

\noindent%
I declare that the work is entirely my own and was produced with no assistance from third parties.\\
I certify that the work has not been submitted in the same or any similar form for assessment to any other examining body and all references, direct and indirect, are indicated as such and have been cited accordingly.
\end{otherlanguage*}


\vspace{3cm}


\noindent%
(\theauthor)\\
Innsbruck, \thedate

\cleardoublepage
}

\pagestyle{fancy}

}
% -----------------------------------------------------------------------------

%
% Self-defined
%
