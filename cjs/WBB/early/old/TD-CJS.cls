\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{TD-CJS}[1994/07/13 v1.2u Standard LaTeX document class]
%
%%%%%%%%%%%%%%%%%%%%% Ignoring Font Warnings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@font@info#1{%
%         \GenericInfo{(Font)\@spaces\@spaces\@spaces\space\space}%
%                     {LaTeX Font Info: \space\space\space#1}
}%
\def\@font@warning#1{%
%         \GenericWarning{(Font)\@spaces\@spaces\@spaces\space\space}%
%                        {LaTeX Font Warning: #1}
}%
\def\@pr@videpackage[#1]{%
  \expandafter\xdef\csname ver@\@currname.\@currext\endcsname{#1}%
  \ifx\@currext\@clsextension
    \typeout{Document Class: \@gtempa\space#1}%
  \else
%    \wlog{Package: \@gtempa\space#1}%
  \fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                DECLARATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newif\ifissue
\global\issuefalse
\DeclareOption{Issue}{\global\issuetrue}



\newif\ifTwelvePT\global\TwelvePTfalse
\DeclareOption{12pt}{\global\TwelvePTtrue}

\newif\ifearlyview
\global\earlyviewfalse
\DeclareOption{earlyview}{\global\earlyviewtrue}

\newif\ifrevise
\global\revisefalse
%
\newif\ifnobookmark
\global\nobookmarkfalse
%
\newif\ifspecial
\global\specialfalse
%
\newif\iftbzdbond\global\tbzdbondfalse
%
\newif\ifmisctitle\global\misctitlefalse %%%
%
\DeclareOption{REVISE}{\global\revisetrue}
\DeclareOption{revise}{\global\revisetrue}
\DeclareOption{MISC}{\global\misctitletrue}
%
\DeclareOption{Thumbnail}{\nobookmarktrue}
%
\newif\ifSDlogo
\global\SDlogofalse
%
\DeclareOption{SDlogo}{\SDlogotrue}
\newif\ifShowThreads
\global\ShowThreadsfalse
\newif\ifHelpNeeded
\global\HelpNeededfalse
%
\DeclareOption{ShowThreads}{%
    \global\ShowThreadstrue
}
\DeclareOption{help}{%
\global\HelpNeededtrue
}
\ProcessOptions
\def\JobName{MODEL 2}
%
\let\DEF\def
\let\DEFINEFont\def
%

\newcount\proofcount
\proofcount0

\newcount\definitioncount
\definitioncount0

\newcount\propositioncount
\propositioncount0

\newcount\corollarycount
\corollarycount0

\newcount\theoremcount
\theoremcount0

\newcount\remarkcount
\remarkcount0

\newcount\notationcount
\notationcount0

\newcount\examplecount
\examplecount0

\newcount\lemmacount
\lemmacount0

\newcount\conjecturecount
\conjecturecount0




\let\NEWCmd\newcommand
\let\RENEWCmd\renewcommand
\let\HANGin\hangindent
\let\GlossaryENV\renewenvironment
%
\newif\ifnameyear
\global\nameyearfalse
%
\newif\iflinenums
\global\linenumsfalse
%
\newif\ifdraft
\global\draftfalse
\newif\iftrim
\newif\ifcontact
\newif\ifprinter
\newif\ifcrop
\global\cropfalse
\newif\ifuncorrectedproof
%\ifcolor\@undefined
\newif\ifcolor
\global\colorfalse
%\else
\global\uncorrectedprooffalse
\newif\ifshowstats
\global\showstatsfalse
\newdimen\draftrule
\newdimen\trimrule
%\draftrule.1pt
%\newdimen\ptdraftrule
%\ptdraftrule0pt
\trimrule.0pt
\newif\ifcountempty
\global\countemptyfalse
\newif\ifnoquery
\global\noqueryfalse
\newif\ifbookreview
\global\bookreviewfalse
\newif\ifrunhead
\global\runheadfalse
\newif\ifhangingindent
\global\hangingindentfalse
\newif\ifIgnoreFigure
\global\IgnoreFigurefalse
\newif\ifequation
\global\equationfalse
\newif\ifeditorial
\global\editorialfalse
\newif\ifnopubdate
\global\nopubdatefalse
%
\newif\ifauto
\global\autofalse
%
%
\newif\ifAbridge\global\Abridgefalse%
\DeclareOption{Abridge}{\global\Abridgetrue}%
%
\newdimen\trimwidth
\newdimen\trimheight
\newdimen\typewidth
\newdimen\typeheight
%
\newdimen\tempdimen
\newdimen\normaltextheight
\newskip\normalbaselineskip
\newskip\belowfrontmatterskip
\belowfrontmatterskip0pt
\newdimen\blindfoliodrop
%
\newdimen\enumdim
\newdimen\figheight
\newdimen\figwidth
\newdimen\tabledim
\newskip\tableleftskip
\newskip\tablerightskip
\newskip\tablebaselineskip
\tablebaselineskip12pt
%
\newdimen\mathindent
\newdimen\emathindent
\newdimen\bibindent
%
\newbox\tempbox
\newbox\tempabox
\newbox\tempbbox
%
\newskip\messageskip
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Language Declarations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\ENGLISH{Language="ENGLISH"}
\def\FRENCH{Language="FRENCH"}
\def\GERMAN{Language="GERMAN"}
\def\SPANISH{Language="SPANISH"}
\def\PORTUGUESE{Language="PORTUGUESE"}
\def\RUSSIAN{Language="RUSSIAN"}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                       COLORS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcount\colorcount
\colorcount1
\def\blue{\special{color cmyk 1 1 0 0}}
\def\white{\special{color cmyk 0 0 0 0.0}}
\def\black{\special{color cmyk 0 0 0 1.0}}
\def\gray{\special{color gray .3}}
\def\yellowtrty{\special{color cmyk 0 .03 .30 0}}
%\def\yellowsxty{\special{color cmyk 0 .06 .60 0}}
\def\red{\special{color cmyk 0 1 \the\colorcount\space 0}}

\long\def\MakeCMYKcolor{\@ifnextchar[{\@CMYKcolorname}{\@@CMYKcolorname}}
\long\def\MakeRGBcolor{\@ifnextchar[{\@RGBcolorname}{\@@RGBcolorname}}

\long\def\@CMYKcolorname[#1]#2{%
\special{color push}\csname CMYK#1\endcsname#2\special{color pop}}

\long\def\@@CMYKcolorname#1#2#3#4#5{%
\special{color push}\special{color cmyk #1 #2 #3 #4}#5\special{color pop}
}

\long\def\@RGBcolorname[#1]#2{%
\special{color push}\csname RGB#1\endcsname#2\special{color pop}}

\long\def\@@RGBcolorname#1#2#3#4{%
\special{color push}\special{color rgb #1 #2 #3}#4\special{color pop}
}

%%%%%%%%%%%%%%%%%%%%%%%% CMYK Color Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\long\def\CMYKCyan{\special{color cmyk 1 0 0 0}}
\long\def\CMYKMeganta{\special{color cmyk 1 0 0 0}}
\long\def\CMYKYellow{\special{color cmyk 1 0 0 0}}
\long\def\CMYKBlack{\special{color cmyk 1 0 0 0}}
\long\def\CMYKRed{\special{color cmyk 0 1 1 0}}
\long\def\CMYKBlue{\special{color cmyk 1 1 0 0}}
\long\def\CMYKGreen{\special{color cmyk 1 0 1 0}}
%
%%%%%%%%%%%%%%%%%%%%%%%% RGB Color Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%\def\RGBRed{\special{color rgb 1 0 0}}
\long\def\RGBRed{\special{color rgb 1 .3 0}}
\long\def\RGBGreen{\special{color rgb 0 1 0}}
\long\def\RGBBlue{\special{color rgb 0 0 1}}
\long\def\RGBBlack{\special{color rgb 0 0 0}}
\long\def\RGBWhite{\special{color rgb 1 1 1}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        OPTIONS DECLARATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\endpage#1{\gdef\@endpage{#1}}\endpage{}%
%
\@twosidetrue\@mparswitchtrue
%
\def\jobtag{Journal Identification = \@jid}
\def\TeXversion{Latex2e}
%
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{nameyear}{\nameyeartrue}
\DeclareOption{BookReview}{\bookreviewtrue}
\DeclareOption{Runhead}{\runheadtrue}
\DeclareOption{Editorial}{\editorialtrue}
\DeclareOption{noquery}{\noquerytrue}
\DeclareOption{countempty}{\countemptytrue}
\DeclareOption{auto}{\global\autotrue}
\DeclareOption{draft}{%
\setlength\overfullrule{0pt}%
\setlength\columnseprule{0.5\p@}%
\setlength\columnseprule{0\p@}%
\draftrule.1pt\drafttrue}
\DeclareOption{showdraft}{%
\setlength\overfullrule{0pt}%
\setlength\columnseprule{0.5\p@}%
\setlength\columnseprule{0\p@}%
\draftrule.1pt}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{contact}{\setlength\overfullrule{0pt}\contacttrue}
\DeclareOption{printer}{\setlength\overfullrule{0pt}\trimrule0pt\printertrue}
\DeclareOption{crop}{\global\croptrue\setlength\overfullrule{0pt}\trimrule0pt\printertrue\draftrule0pt\draftfalse}
\DeclareOption{linenums}{\linenumstrue}
\DeclareOption{showtrimmarks}{%
%  \@draftfalse
  \def\query{}%
  \overfullrule \z@
\trimrule.1pt\trimtrue
}
\DeclareOption{hidetrimmarks}{%
%  \@draftfalse
  \def\query{}%
  \overfullrule \z@
\trimrule0pt\trimtrue
}
\DeclareOption{uncorrectedproof}{%
\uncorrectedprooftrue
}
\DeclareOption{showstats}{%
\showstatstrue
}
\DeclareOption{ignorefigures}{%
\IgnoreFiguretrue
}
\DeclareOption{ignfigs}{%
\IgnoreFiguretrue
}
\DeclareOption{colorlinks=true}{%
  \colortrue
}
\DeclareOption{colorlinks=false}{%
  \colorfalse
}

\DeclareOption{special}{%
  \specialtrue
}
%
\newif\ifArticleInPress
\global\ArticleInPressfalse

\newif\ifContinuous
\global\Continuousfalse
%
\newif\ifSecContinuous
\global\SecContinuousfalse
%
\newcounter{ContinuousCnt}%
\setcounter{ContinuousCnt}{0}%
\renewcommand\theContinuousCnt{\arabic{ContinuousCnt}}
%
\newif\ifSeqEnun
\global\SeqEnunfalse
%
\newif\ifSecSeq
\global\SecSeqfalse


\DeclareOption{Enun="Cont"}{%
    \global\Continuoustrue}

\DeclareOption{Enun="Seq"}{%
          \global\SeqEnuntrue}

\DeclareOption{Enun="SecSeq"}{%
\global\Continuousfalse%
\global\SecContinuousfalse%
\global\SeqEnunfalse%
\global\SecSeqtrue
}%


\DeclareOption{Enun="SecCont"}{%
\renewcommand\theContinuousCnt{\thesection.\arabic{ContinuousCnt}}
\@addtoreset{ContinuousCnt}{section}
\Continuoustrue\global\SecContinuoustrue
}

%

\def\CropPg#1#2#3#4{\gdef\@CropPg{#1 #2 #3 #4}}
\CropPg{162.5}{237}{647}{957}
%
\PassOptionsToPackage{dviwindo}{rotating}
\ExecuteOptions{final}%,ignfigs
\ProcessOptions
%
\usepackage{dcolumn}
%\newcolumntype{.}{D{.}{.}{10}}
\newcolumntype{,}{D{,}{,}{-1}}
\newcolumntype{(}{D{(}{}{-1}}
\newcolumntype{.}{D{.}{.}{1}}
%
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\bibname{REFERENCES}
\newcommand\indexname{Index}
\newcommand\figurename{Fig}
\newcommand\tablename{\sc{Table}}
\newcommand\partname{Part}
\newcommand\chaptername{Chapter}
\newcommand\appendixname{Appendix}
%
\newcommand\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
%
\newcount\minute
\newcount\hour
\def\currenttime{%
     \minute\time
     \hour\minute
     \divide\hour60
\ifnum\hour=12
     12:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=13
     1:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=14
     2:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=15
     3:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=16
     4:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=17
     5:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=18
     6:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=19
     7:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=20
     8:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=21
     9:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=22
     10:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=23
     11:\multiply\hour60\advance\minute-\hour\the\minute\,pm
\else
\ifnum\hour=0
     12:\multiply\hour60\advance\minute-\hour\the\minute\,am
\else
     \the\hour:\multiply\hour60\advance\minute-\hour\the\minute\,am
\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           FONTS DECLARATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\DEFINEFont\rhfont{\sffamily\fontsize{9}{11}\selectfont}%
\DEFINEFont\foliofont{\sffamily\fontsize{9}{11}\selectfont}%
\DEFINEFont\dropfoliofont{\fontsize{9}{9}\selectfont}%


%
\DEFINEFont\PNfont{}%
\DEFINEFont\PTfont{}%
%
\DEFINEFont\reprintfont{\fontsize{7}{10}\selectfont}
\DEFINEFont\historyfont{\fontsize{10}{12}\selectfont}
\DEFINEFont\journalurlfont{\fontsize{8}{10}\selectfont}
\DEFINEFont\documentheadfont{\fontsize{17}{18}\selectfont\leftskip0pc \rightskip0pc plus1fill}
\DEFINEFont\authorfont{\sffamily\ifTwelvePT\fontsize{11}{13}\selectfont\else\fontsize{10}{12}\selectfont\fi\leftskip0pc \rightskip0pc\let\it\itshape}
\DEFINEFont\affiliationfont{\ifTwelvePT\fontsize{10}{13}\selectfont\else\fontsize{9}{12}\selectfont\fi\itshape\raggedright}
\DEFINEFont\CNfont{\fontsize{17}{20}\selectfont\bfseries}%\mathversion{bold}
\DEFINEFont\CTfont{\sffamily\fontsize{17}{19}\selectfont\bfseries\leftskip0pc \rightskip1pc plus1fill\let\it\itshape}%
\DEFINEFont\ECTfont{\sffamily\fontsize{17}{19}\selectfont\bfseries\raggedright}%
\DEFINEFont\subtitlefont{\fontsize{17}{20}\selectfont\leftskip0pc plus1fill\rightskip0pc plus1fill}
\DEFINEFont\summaryfont{\ifTwelvePT\fontsize{10}{12}\selectfont\else\fontsize{9}{11}\selectfont\fi\let\it\itshape\let\bf\bfseries}
\DEFINEFont\abstractfont{\ifTwelvePT\fontsize{10}{12}\selectfont\else\fontsize{9}{11}\selectfont\fi\let\it\itshape\let\bf\bfseries}
\DEFINEFont\summarysectionfont{\sffamily\fontsize{10}{12}\selectfont}
\DEFINEFont\Keyheadfont{\sffamily\fontsize{10}{12}\selectfont}
\DEFINEFont\Keywordsfont{\fontsize{10}{12}\selectfont}
%
\DEFINEFont\footcopyfont{\fontsize{8}{10}\selectfont}%
\DEFINEFont\piifont{\fontsize{8}{8}\selectfont}%
%
\DEFINEFont\sectionfont{\let\it\itshape\sffamily\fontsize{10}{12}\selectfont\rightskip1pc plus 1fill\let\rmbeta\brmbeta\let\rmalpha\brmalpha\let\rmgamma\brmgamma\let\rmdelta\brmdelta\let\rmtheta\brmtheta\let\rmmu\brmmu\let\rmkappa\brmkappa\let\rmeta\brmeta\ifappendix\global\addtocounter{appcount}{1}\fi}%\mathversion{bold}
\DEFINEFont\subsectionfont{\let\bf\bfseries\sffamily\fontsize{10}{12}\selectfont\let\it\itshape\rightskip1pc plus 1fill}%
\DEFINEFont\subsubsectionfont{\let\bf\bfseries\sffamily\sffamily\fontsize{10}{12}\selectfont\itshape\let\it\itshape\let\sf\SubSecsf\rightskip1pc plus 1fill\let\rmPsi\itPsi\let\rmlambda\itlambda\let\rmsigma\itsigma\let\rmbeta\itbeta\let\rmalpha\italpha\let\rmgamma\itgamma\let\rmdelta\itdelta\let\rmtheta\ittheta\let\rmmu\itmu\let\rmphi\itphi\let\Gamma\itGamma\let\Delta\itDelta\let\Theta\itTheta\let\Lambda\itLambda\let\Phi\itPhi\let\Psi\itPsi\let\Omega\itOmega\let\Xi\itXi\let\Pi\itPi\let\Sigma\itSigma\let\Upsilon\itUpsilon\let\rmkappa\itkappa\let\rmpi\itpi\let\rmeta\iteta}%%
\DEFINEFont\paragraphfont{\sffamily\fontsize{10}{12}\selectfont\bfseries\let\bf\bfseries\let\it\itshape}%
\DEFINEFont\subparagraphfont{\subsectionfont}%
\DEFINEFont\subsubparagraphfont{}%
\DEFINEFont\MandMsectionfont{\fontsize{9}{11}\selectfont\bfseries\rightskip1pc plus 1fill}%\mathversion{bold}
\DEFINEFont\MandMsubsectionfont{\fontsize{9}{11}\selectfont\itshape\rightskip1pc plus 1fill}%
%
\DEFINEFont\figcaptionfont{\fontsize{9}{10}\selectfont\let\it\itshape}
\DEFINEFont\figcaptionnumfont{\fontsize{9}{10}\selectfont}%\bfseries\scshape
%
\DEFINEFont\tablefont{\fontsize{9}{12}\selectfont}%\leftskip\tableleftskip\rightskip\tablerightskip}% use plus 1fill if needed
\DEFINEFont\tablecaptionfont{\fontsize{9}{12}\selectfont\leftskip0pt plus1fill\rightskip0pt plus1fill}%
\DEFINEFont\tablecaptionnumfont{\fontsize{9}{11}\selectfont}%
\DEFINEFont\TCHfont{}%
\DEFINEFont\TSCHfont{}%
\DEFINEFont\tabnotefont{\fontsize{9}{11}\selectfont\leftskip0pt plus1fill\rightskip0pt plus1fill}% use plus 1fill if needed
%
\DEFINEFont\indexfont{\footnotesize\raggedright}
\DEFINEFont\bibliofont{\ifTwelvePT\fontsize{10}{12}\selectfont\else\fontsize{9}{11}\selectfont\fi}
%
\DEFINEFont\listfont{}
\DEFINEFont\listdevicefont{}
\DEFINEFont\extractfont{\leftskip1pc}
\DEFINEFont\quotefont{\leftskip2pc}
\DEFINEFont\exhheadfont{\fontsize{9}{11}\selectfont\bfseries}
\DEFINEFont\exhhead#1{{\exhheadfont #1\par}\@afterheading}
\DEFINEFont\leadintextfont{\itshape}
\DEFINEFont\leadintext#1{{\leadintextfont#1}\hskip2.5pt}%
\DEFINEFont\lhfont{\fontsize{10}{12}\selectfont\bfseries}
\DEFINEFont\lh#1{\par\addvspace{\baselineskip}\noindent{\lhfont#1\par}\vspace*{-5.9pt}}
\DEFINEFont\verbatim@Font{\itshape\normalfont}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

\def\Arno#1{\gdef\@Arno{#1}}
\Arno{\@jid\ \@aid}

\def\endmark{}
%
\def\spreadoutfactor{0}
%
\def\spreadout#1{%
     \gdef\temp{#1}%
  \dimen0 = \spreadoutfactor pt
  \ifdim\dimen0=\z@\temp\else                            % Modified by Tarun. If the spreadoutfactor
  \expandafter\dospreadout\temp\endmark\kern-\dimen0\fi} % is zero, then we can use kern or hskip
                                                         % with the title etc. Suppressed the extra space
                                                         % after spreadout.
\def\dospreadout{%
  \afterassignment\findospreadout
  \let\next= }
%
\def\findospreadout{%
   \ifx\next\endmark
      \let\nextaction = \relax
   \else
      \let\nextaction = \dospreadout
      \next
      \kern\dimen0
   \fi
   \nextaction}
%
\def\@xvpt{15}
\def\@xviiipt{18}
\def\@xxivpt{24}
%
\normalbaselineskip12pt
%
\newif\ifMaterialsandMethods
\global\MaterialsandMethodsfalse
%
\renewcommand\normalsize{%
\ifTwelvePT
   \@setfontsize\normalsize{12}{14}%\@xipt{13}
\else
   \@setfontsize\normalsize{10}{12}%\@xipt{13}
\fi
   \abovedisplayskip 10\p@ \@plus2\p@ \@minus.5\p@
   \abovedisplayshortskip \abovedisplayskip
   \belowdisplayshortskip \abovedisplayskip
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\def\biggg#1{{\hbox{$\left#1\vbox to20.5\p@{}\right.\n@space$}}}
\def\Biggg#1{{\hbox{$\left#1\vbox to23.5\p@{}\right.\n@space$}}}
\normalsize
%
\newcommand\small{%
   \@setfontsize\small{8}{10}%\@ixpt{11}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setSmallDelims}
%
\def\setSmallDelims{%
\def\big##1{{\hbox{$\left##1\vbox to7.5\p@{}\right.\n@space$}}}%
\def\Big##1{{\hbox{$\left##1\vbox to10.5\p@{}\right.\n@space$}}}%
\def\bigg##1{{\hbox{$\left##1\vbox to13.5\p@{}\right.\n@space$}}}%
\def\Bigg##1{{\hbox{$\left##1\vbox to16.5\p@{}\right.\n@space$}}}%
\def\biggg##1{{\hbox{$\left##1\vbox to19.5\p@{}\right.\n@space$}}}%
\def\Biggg##1{{\hbox{$\left##1\vbox to22.5\p@{}\right.\n@space$}}}%
}
%
\newcommand\footnotesize{%
   \@setfontsize\footnotesize{9}{10}%\@viiipt{10}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setFootnotesizeDelims}
%
\newcommand\FFootnotesize{%
   \@setfontsize\footnotesize{10}{12}%\@viiipt{10}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setFootnotesizeDelims}
%
\def\setFootnotesizeDelims{%
\def\big##1{{\hbox{$\left##1\vbox to6.5\p@{}\right.\n@space$}}}%
\def\Big##1{{\hbox{$\left##1\vbox to9.5\p@{}\right.\n@space$}}}%
\def\bigg##1{{\hbox{$\left##1\vbox to12.5\p@{}\right.\n@space$}}}%
\def\Bigg##1{{\hbox{$\left##1\vbox to15.5\p@{}\right.\n@space$}}}%
\def\biggg##1{{\hbox{$\left##1\vbox to18.5\p@{}\right.\n@space$}}}%
\def\Biggg##1{{\hbox{$\left##1\vbox to21.5\p@{}\right.\n@space$}}}%
}
%
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% For Times family %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareOldFontCommand{\bi}{\bfseries\itshape}{\bfseries\itshape}
\newcommand{\cal}{\protect\pcal}
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
\newcommand{\mit}{\protect\pmit}
\newcommand{\pmit}{\@fontswitch{\relax}{\mathnormal}}





%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\trimheight{254mm}%
\setlength\trimwidth{171mm}%
%
\setlength\typeheight{50pc}%
\setlength\typewidth{135mm}%
%
\newlength\extramargin
\setlength\extramargin{0pc}
%
\setlength\textheight{50pc}%
\addtolength\textheight{-4.85pt}
\setlength\textwidth{135mm}
\setlength\normaltextheight{50pc}
\addtolength\normaltextheight{-4.85pt}
%
\setlength\topmargin{13.2mm}%
\setlength\oddsidemargin{18mm}%
\setlength\evensidemargin\trimwidth
\addtolength\evensidemargin{-\textwidth}
\addtolength\evensidemargin{-\oddsidemargin}
%
\setlength\headheight{7\p@}% 6.4
\setlength\headsep   {17pt}% \typeheight - \textheight - \headheight
\setlength\topskip   {0\p@}
%
\setlength\footskip{11.3pt}
\setlength\maxdepth{4\p@}
%
\blindfoliodrop\trimheight
\advance\blindfoliodrop-\typeheight
\advance\blindfoliodrop-\topmargin
\advance\blindfoliodrop-\footskip
\advance\blindfoliodrop18pt
%
\setlength\parindent{15pt}%
%
\setlength\marginparwidth {.75in}
\setlength\marginparsep{7\p@}
\setlength\marginparpush{5\p@}
%
\setlength\footnotesep{9\p@}%6.65
\setlength{\skip\footins}{18\p@ \@plus 4\p@ \@minus 3\p@}% change + based on FN font height
\skip\@mpfootins = \skip\footins
%
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 0\p@}
\setlength\textfloatsep{11\p@ \@plus 4\p@ \@minus 0\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 0\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 0\p@}
\setlength\dbltextfloatsep{15\p@ \@plus 4\p@ \@minus 0\p@}

% FOR FLOATS ON A SEPARATE FLOAT PAGE OR COLUMN:
% ONE-COLUMN MODE OR SINGLE-COLUMN FLOATS IN TWO-COLUMN MODE:
\@fptop 0pt plus 1fil    % Stretch at top of float page/column. (Must be
                         % 0pt plus ...)
\@fpsep 8pt plus 2fil    % Space between floats on float page/column.
\@fpbot 0pt plus 1fil    % Stretch at bottom of float page/column. (Must be
                         % 0pt plus ... )
%
% DOUBLE-COLUMN FLOATS IN TWO-COLUMN MODE.
\@dblfptop 0pt plus 1fil % Stretch at top of float page. (Must be 0pt plus ...)
\@dblfpsep 8pt plus 2fil % Space between floats on float page.
\@dblfpbot 0pt plus 1fil % Stretch at bottom of float page. (Must be
                         % 0pt plus ... )
%
\setlength\partopsep{0pt}
\setlength\lineskip{1\p@}% check if it can be flexible
\setlength\normallineskip{1\p@}%
\renewcommand\baselinestretch{}
\setlength\parskip{1pt}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
%
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty=1000%     -\@lowpenalty
\clubpenalty=10000
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                PART TITLE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand\part{\cleardoublepage
                    \thispagestyle{empty}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@part\@spart}
%
\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    \markboth{}{}
     \ifnum \c@secnumdepth >-2\relax
       \noindent{\PNfont \partname~\thepart\par}%
       \vskip 20\p@
     \fi
     \noindent{\PTfont #2\par}%
    \@endpart}
\def\@spart#1{%
    {\centering
     \interlinepenalty \@M
     \Huge \bfseries #1\par}%
    \@endpart}
%
\def\@endpart{\vfil\newpage}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                FRONT MATTER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%\input epsf.def
\usepackage{epsfig}

\def\@ins@ps#1{\special{psfile="#1"}}
\def\subtitle#1{%
    \gdef\@subtitle{#1}}
\subtitle{}%
\def\jvol#1{\gdef\@jvol{{#1}}}
\jvol{}
\def\issue#1{\gdef\@issue{#1}}
\issue{}


\def\jid#1{\gdef\@jid{#1}}%
\jid{}
\def\rhjvol#1{\gdef\@rhjvol{{{#1}}}}\rhjvol{}
\def\rhissue#1{\gdef\@rhissue{{{#1}}}}\rhissue{}
\def\aid#1{\gdef\@aid{#1}}\aid{}
\def\jyear#1{\gdef\@jyear{#1}}
\def\jname#1{\gdef\@jname{\it #1}}
\jname{\csname\@jid\endcsname}


\def\cjyear#1{\gdef\@cjyear{#1}}
\cjyear{}

\def\Fjname#1{\gdef\@Fjname{\it #1}}
\Fjname{La revue canadienne de statistique}

\jyear{2020}
\jvol{}
\jname{The Canadian Journal of Statistics}
\def\PrintReprintLine{}

\def\Fjname#1{\gdef\@Fjname{\it #1}}
\Fjname{La revue canadienne de statistique}



\def\rhauthor#1{\gdef\@rhauthor{\uppercase{#1}}}\rhauthor{}
\def\rhtitle#1{\gdef\@rhtitle{\uppercase{#1}}}\rhtitle{}


\def\Twolinereprint#1{%
    \gdef\@Twolinereprint{\centering #1\\ \ifx\@jvol\@empty\@aid\else\@jvol\fi\ \@jyear\ }}
\Twolinereprint{\@jname}

\def\pubdate#1{{\piifont\gdef\@pubdate{#1}}\par}%
\pubdate{}

\def\received#1{%
\gdef\@received{#1}}
\received{}%


\def\revised#1{%
\gdef\@revised{#1}}
\revised{}%

\def\accepted#1{%
\gdef\@accepted{#1}}
\accepted{}


\def\journalurl#1{\gdef\@journalurl{#1}}

\journalurl{www.elsevier.nl/locate/optmat}

\def\documenthead#1{\gdef\@documenthead{#1}}\documenthead{}
\def\fnm#1{#1}
\def\snm#1{{#1}}
\def\degs#1{#1}
\let\dochead\documenthead
\let\atlhead\documenthead

\def\issnnum#1{\gdef\@issnnum{#1}}\issnnum{}





\def\copyrightline#1{%
\newif\ifabscop
\global\abscopfalse
\gdef\@copyrightline{\ifabscop\else\par\fi#1}
\gdef\@checkcopyrightline{#1}
\gdef\@copyrightabs{#1}}
\copyrightline{}%will change in abstract also


\def\Frenchcopyrightline#1{%
\newif\ifabscop
\global\abscopfalse
\gdef\@Frenchcopyrightline{\ifabscop\else\par\fi#1}
\gdef\@checkFrenchcopyrightline{#1}
\gdef\@Frenchcopyrightabs{#1}}
\Frenchcopyrightline{}


\def\PIInumber#1{\gdef\@forPDFTitle{PII: #1}%
\gdef\@PIInumber{PII: \def\spreadoutfactor{1}\spreadout{#1}}}\PIInumber{}

\def\DOInumber#1{\gdef\@forPDFTitle{DOI: #1}%
\gdef\@DOInumber{\par\hfill{\piifont DOI: #1\par}}}\DOInumber{}%

\def\copyrightStatus#1{\gdef\@copyrightStatus{#1}}
\copyrightStatus{}

\gdef\@STAT{004}
\gdef\@STATA{003}


\def\@copyright{\vbox{\hbox{\footcopyfont%
{\fontsize{9}{9}\selectfont\sf{$\Ccopyright$\ \@jyear\ John Wiley \& Sons, Ltd.}}%
}}}

\def\@catchline{%
\lower0pt\vbox{\reprintfont\hbox{{\@jname}\par}
\ifearlyview\else\hbox{Vol. \@jvol\ifx\@issue\@empty\else,\ No. \@issue\fi,\ \@jyear, Pages\ \ifx\@firstpage\@lastpage{\@firstpage}\else\@firstpage--\@lastpage\fi}\par\fi
\hbox{\@Fjname}\par}}

\def\abshead#1{\vskip8pt\par\noindent{\fontsize{9}{11}\selectfont\itshape #1}}

\newif\ifnoabs
\global\noabsfalse

\def\keywords#1{{\par\vspace*{6.4pt}
\ifTwelvePT\fontsize{10}{12}\selectfont\else\fontsize{9}{11}\selectfont\fi\let\qtoa\Qtoa%
\def\KWD{{\vskip9.5pt\par\it Keywords: }}
\def\KWDtitle##1{{\vskip3.9pt\par\parindent0pt\it ##1:\hskip3pt}}
\def\FKWDtitle##1{{\vskip7.5pt\par\parindent0pt\it ##1:\hskip3pt}}
\def\MSC{\vskip6pt\par\leavevmode\hbox {\it MSi:\ }}%
\def\AMS{\vskip6pt\par\leavevmode\hbox {\it AMS Classification:\ }}%
\def\PACS{\vskip6pt\par\leavevmode\hbox {\it PACS:\ }}%
\def\IDT##1{\vskip6pt\par\leavevmode\hbox {##1:\ }}%
\def\JEL{\vskip6pt\par\leavevmode\hbox {\it JEL classification:\ }}%
\def\UK{\vskip6pt\par\leavevmode\hbox {\textit{Keywords:\ }}}%
\def\ABR{{\vskip6pt\par\it Abbreviations: }}
\let\CK\relax
#1\par\addvspace{1.3pt}}}
\def\sep{;\ }
\def\fsep{ ;\ }




\newcount\myabseccount
\myabseccount0
\gdef\abstractsection#1{%
\gdef\@checkabssec{#1}
\ifnum\myabseccount>0
    \par\noindent\hbox{\itshape #1:\hskip3pt\ignorespaces}%
    \else%
    {\ifx\@checkabssec\@empty\else\noindent #1\fi\par\vskip6.5pt}%
\fi\advance\myabseccount1}

\newcount\myfabseccount
\myfabseccount0
\gdef\fabstractsection#1{%
\gdef\@checkfabsec{#1}%
\ifnum\myfabseccount>0
    \par\noindent\hbox{\itshape #1:\hskip3pt\ignorespaces}%
    \else%
    {\par\ifx\@checkfabsec\@empty\else\addvspace{11.5pt}\fi\noindent \fontsize{10}{12}\selectfont #1\par\vskip2.7pt}%
\fi\advance\myfabseccount1}

\def\abstractsubsection#1{{\vskip11pt\itshape\noindent #1.}}

\def\Frenchcopyright{\par%
\noindent\hbox{\@copyrightline}\global\myfabseccount0\global\myabseccount0\par\vspace*{6pt}}

\newif\ifAbstract
\global\Abstractfalse

\newif\ifabszdbond
\global\abszdbondfalse

\newif\ifIndent
\global\Indentfalse
\def\Abstract{\global\Abstracttrue\global\Indenttrue}
\def\Indent{\global\Indenttrue}

\def\abstract{\@ifnextchar[{\@abstract}{\@abstract[]}}
\def\@abstract[#1]{\global\Abstracttrue\abszdbondtrue%
\def\LANGUAGE{#1}
\setlength\parindent{9pt}%
%    \noindent\vrule width\textwidth height.25pt depth0pt
    \vskip5pt
        \abstractfont
    \ifx\LANGUAGE\ENGLISH
            \noindent{\bfseries Abstract}\par
        \vskip5pt
    \else
    \ifx\LANGUAGE\FRENCH
            \noindent{\bfseries Résumé}\par
        \vskip5pt
    \else
    \ifx\LANGUAGE\GERMAN
            \noindent{\bfseries Zusammenfassung}\par
        \vskip5pt
    \else
    \ifx\LANGUAGE\SPANISH
            \noindent{\bfseries Resumen}\par
        \vskip5pt
    \else
    \ifx\LANGUAGE\PORTUGUESE
            \noindent{\bfseries Sumário}\par
        \vskip5pt
    \else
    \ifx\LANGUAGE\RUSSIAN
            \noindent{\bfseries peepat}\par
        \vskip5pt
    \fi\fi\fi\fi\fi\fi}
%
\def\endabstract{\par%
\let\Qtoa\qtoa}
%
\long\def\startabstract#1{%
    \gdef\@startsummary{#1%
        \par
\ifAbstract
        \vskip0.2pt%
\else
    \vskip1pt
\fi
\vskip\belowfrontmatterskip
\fi}}
\startabstract{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Author & Affiliation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcount\authorcount
\authorcount0
\newcount\affiliationcount
\affiliationcount0
\newcount\tempcount
\newcount\affiltempcount
\newcount\checkauthorcount
\checkauthorcount=0
\newcount\authorendcount

\newcount\plncount
\plncount0
\newcount\plntempcount
%
\newif\ifbreak
\global\breakfalse
%
%

%%%%%%%%%%%%%%%%%%%%%
%
\newcount\acount
\acount=0
\newcount\auttemp

\IfFileExists{\jobname.aut}{\input{\jobname.aut}}{%
%\message{^^J%
%*********************************************************^^J%
%Recompile the file for generation of author names.^^J%
%*********************************************************^^J%
%^^J}
%\read-1 to \messageskip
}%
%
%%%%%%%%%%%%%%%%%%%%%
%
\def\author#1{%
    \global\advance\authorcount by 1%
    \immediate\write\aut{\string\auttemp=\the\authorcount}%
    \acount\auttemp
    \advance\acount-2
        \ifnum\authorcount=\acount
            \expandafter\gdef\csname author\romannumeral\the\authorcount\endcsname{#1\vphantom{p} and}%
        \else
        \ifnum\authorcount>\acount
            \expandafter\gdef\csname author\romannumeral\the\authorcount\endcsname{#1\vphantom{p}}%
        \else
            \expandafter\gdef\csname author\romannumeral\the\authorcount\endcsname{#1\vphantom{p},}%
        \fi\fi%
}
\newwrite\aut
\immediate\openout\aut=\jobname.aut




\def\corfont{\sffamily\fontsize{10}{12}\selectfont}

\def\docheadfootnote#1{{%\renewcommand\thefootnote{}%
\insert\footins{\footnotesize\leavevmode\hbox to 15.8pt{\hfill$^{1}$\hskip2.5pt}#1}}}%




\def\email#1{%
{{\it E-mail: #1}}}


\long\def\@MakeFnText#1{\overfullrule0pt%
\noindent\smash{\@makefnmark}{#1}}

\def\correspondingauthor[#1]#2{%
\footnotetext{\vbox{\hsize\textwidth{\hspace*{-1.7pt}\footnotesize{{#1\,}\hypertarget{cor1}{\it Author to whom correspondence may be addressed.}}\ #2}}\setcounter{footnote}{1}}}%
%\insert\footins{\hsize\textwidth{\hspace*{-10pt}\footnotesize{{#1\,}\hypertarget{cor1}{\it Author to whom correspondence may be addressed.}}\ #2\par}\vspace*{-1.5pt}}\setcounter{footnote}{1}}%

\newcount\contractcount
\contractcount0



\def\contract#1{%
\insert\footins{%
\global\advance\contractcount1
\hsize\textwidth{\ifnum\contractcount>1\vskip1pt\else\vskip10pt\fi\parindent0pt\footnotesize{Contract/grant sponsor:}\ #1\par}}}%


\def\contracts#1{%
\insert\footins{%
\global\advance\contractcount1
\hsize\textwidth{\ifnum\contractcount>1\vskip1.8pt\else\vskip10pt\fi\parindent0pt\footnotesize{Contract/grant sponsors:}\ #1\par}}}%

\def\sbreak{\vspace*{-1.8pt}\break}
\def\sibreak{\vspace*{-1.4pt}\break}
\def\siibreak{\vspace*{-1.2pt}\break}
\def\siiibreak{\vspace*{-1.8pt}\break}
\def\tbreak{\vspace*{-4.2pt}\break}
\def\abreak{\vspace*{-1.6pt}\break}

\newcount\Authcount
\Authcount0

\def\authorref#1{%
\advance\Authcount1
\ifnum\Authcount>1
\smash{$^{\hyperlink{affil#1}{{\hbox{\scriptsize{#1}}}}}$}%
\else
\smash{$^{\hyperlink{affil#1}{{\hbox{\scriptsize{#1}}}}}$}%
\fi}


\newcount\thankscount
\thankscount0

\def\thanksref#1{%
\advance\thankscount1
\ifnum\thankscount>1
\smash{${\hyperlink{cor1}{\hbox{{#1}}}}$}%
\else
\smash{${\hyperlink{cor1}{\hbox{{#1}}}}$}%
\fi}

\newcounter{authorfootnotecounter}
\setcounter{authorfootnotecounter}{0}
\def\authorfootnote#1#2{\renewcommand\thefootnote{}%
\stepcounter{authorfootnotecounter}
\footnotetext[-1]{\hspace*{-6pt}$^{#1}$\hspace*{3.5pt}#2}}
%
\def\affiliation{\@ifnextchar[{\@@affiliation}{\@@@affiliation}}
\def\@@affiliation[#1]#2{\global\advance\affiliationcount by 1
   \expandafter\gdef\csname affiliation\romannumeral\the\affiliationcount\endcsname{$^{{\hbox{\scriptsize\selectfont{\rm{#1}\hypertarget{affil#1}{}}}}}$#2\vphantom{y}\par}}
\def\@@@affiliation#1{\global\advance\affiliationcount by 1
   \expandafter\gdef\csname affiliation\romannumeral\the\affiliationcount\endcsname{{}#1\par}}
%
\def\Hrulefill{\leavevmode\leaders\hrule \@height4.5pt depth -4.25pt\hfill}

%\newcount\cskip

\def\makeauthor{%\setlength\parindent{0pt}%
\def\snm##1{{##1}}
   \tempcount0%
   \vbox{\hsize\textwidth\authorfont\raggedright
   \loop%
      \ifnum\tempcount<\authorcount%
         \advance\tempcount by 1%
         \csname author\romannumeral\the\tempcount\endcsname\ %
   \repeat%
\par}\vskip-5pt}
%
\def\makeaffiliation{\setlength\parindent{10pt}%
   \tempcount0%
   {\affiliationfont%
   \loop%
      \ifnum\affiltempcount<\affiliationcount%
         \advance\affiltempcount by 1%
         \csname affiliation\romannumeral\the\affiltempcount\endcsname%\endgraf
 \repeat
\par}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                CHAPTER TITLE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\cleardoublepage{\clearpage\if@twoside\ifodd\c@page\else
    \hbox{}\thispagestyle{empty}\newpage\if@twocolumn\hbox{}\newpage\fi\fi\fi}
%
\newcommand\@chapapp{}%\chaptername}
%
\newcommand\chapter{\clearpage%\cleardoublepage
                    \thispagestyle{plain}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@chapter\@schapter}
%
\def\@chapter[#1]#2{\refstepcounter{chapter}%
                    \gdef\@chaptermark{\uppercase{#1}}%
                    \gdef\@chaptertitle{{{#2}}}}
\def\reviewname{Book Review}
\def\reviewfont{\fontsize{12}{12}\selectfont\leftskip0pt plus1fill\rightskip\leftskip}

\def\collabauthor#1{\gdef\@collabauthor{#1}}
\collabauthor{}
\def\collabaffiliation#1{\gdef\@collabaffiliation{#1}}
\collabaffiliation{}





\newbox\Mybox

\newif\ifRefersdoi
\newif\ifReferspii
\newif\ifRefersDois %%%% 25-02-05
\newif\ifRefersPiis


\def\Refersdoi{\Referspiifalse\Refersdoitrue}
\def\Referspii{\Referspiitrue}

\def\RefersDois{\RefersPiisfalse\RefersDoistrue}
\def\RefersPiis{\RefersPiistrue}


\def\refersdoi#1{\Refersdoitrue%
\gdef\@refersdoi{#1}}

\def\referspii#1{\Referspiitrue%
\gdef\@referspii{#1}}


\def\RefersDois#1{\RefersDoistrue%
\gdef\@RefersDois{#1}}

\def\RefersPiis#1{\RefersPiistrue%
\gdef\@RefersPiis{#1}}


\def\Spchar#1{\gdef\@Spchar{\foliofont#1}}
\Spchar{S}


\def\@makechapterhead#1{%\leftskip0pt plus1fill\rightskip\leftskip%
\setlength\parskip{0pt}
\setlength\parindent{0pt}
\ifnum \c@secnumdepth >\m@ne%
\vbox to 0pt{\vskip25.5pt\documentheadfont{\@documenthead}}%\vskip0\baselineskip
\ifx\@documenthead\empty\vskip19.5pt\else\vskip50.7pt\fi
\ifeditorial
\vspace*{0.5pt}
\noindent{\ECTfont{}#1\vphantom{y}\par\vspace*{33.5pt}}\else\fi%
\ifeditorial\else\noindent{\CTfont{}{\textsf#1}\vphantom{y}\ifx\@subtitle\@empty\else:\quad\@subtitle\fi\par}\fi%
\ifnum\authorcount<1\vskip-9.7pt\else\vskip12.8pt\fi
\makeauthor
\ifnum\affiliationcount<1\vspace*{5pt}\else\vspace*{15pt}\fi
\makeaffiliation
\vskip.5pt
\@startsummary\setlength\parskip{1pt}}

\def\@schapter#1{\chaptermark{#1}%
                 \if@twocolumn
                   \@topnewpage[\@makeschapterhead{#1}]%
                 \else
                   \@makeschapterhead{#1}%
                   \@afterheading
                 \fi}

\def\@makeschapterhead#1{%
  \vbox{\raggedright%\parindent \z@
    \noindent{\CTfont #1}\par
    \nobreak
    \vskip 20\p@
  }}
%
\let\title\chapter
%
\countdef\L@page=0 \L@page=1
\def\makechaptertitle{\thispagestyle{copyright}%
                    \if@twocolumn
                   \twocolumn[%
                   \@makechapterhead{\@chaptertitle}%
\vspace*{2.3pc}%
                   \vskip\belowfrontmatterskip%
]%
                   \else
                   \@makechapterhead{\@chaptertitle}%
                   \fi
                   \gdef\@author{}%
                   %\ifIndent\else\@afterheading\fi%\@afterheading%
                    \setlength\parindent{15pt}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Chapter Title %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Page Numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%\AtBeginDocument{\input{\jobname.pag}}
\def\firstpage#1{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \gdef\@firstpage{#1}%\gdef\@lastpage{#1}%
  \global\c@page=#1\ignorespaces\fi
  }
\def\@firstpage{1}

\def\lastpage#1{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \gdef\@lastpage{#1}\ignorespaces\fi}

\def\@lastpage{1}
\def\@pagerange{1--1}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                SECTION COMMANDS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{secnumdepth}{4}
\newcounter {part}
\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thepart          {\Roman{part}}
\renewcommand\thechapter       {\arabic{chapter}}
\renewcommand\thesection       {\arabic{section}}%\thechapter.
\renewcommand\thesubsection    {\thesection .\arabic{subsection}}%{\thesection.\Alph{subsection}}
\renewcommand\thesubsubsection {\thesubsection .\arabic{subsubsection}}
\renewcommand\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\renewcommand\thesubparagraph  {\theparagraph.\arabic{subparagraph}}
%
%%%%%%%%%%%%% For the generation of section Number %%%%%%%%%%%
%
\def\Secno#1{#1.\hskip5pt\ignorespaces}%{\xdef\@Secno{#1}}\Secno{}
%
%\def\@seccntformat#1{\csname the#1\endcsname.\hskip4.3pt}%\hskip 1em
\def\@subseccntformat#1{\csname the#1\endcsname.\hskip3.3pt}%\hskip 1em
\def\@prgformat#1{\csname the#1\endcsname.\hskip0pt}%\hskip 1em
%
\newcommand\section{\@startsection {section}{1}{\z@}{-17pt plus-1pt minus-1pt}{3pt}{\sectionfont}}%
\newcommand\subsection{\@startsection{subsection}{2}{\z@}{-8.7pt plus-3pt minus-.2pt}{1.9pt}{\subsectionfont}}%
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{-8.7pt plus-3pt minus-.2pt}{1.5pt}{\subsubsectionfont}}%
\newcommand\paragraph{\@startsection{paragraph}{4}{\parindent}{0pt plus3pt minus.2pt}{-9.5pt}{\paragraphfont}}%
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}{-12pt plus-3pt minus-.2pt}{-6pt}{\subparagraphfont}}%
\newcommand\subsubparagraph{\@startsection{subsubparagraph}{6}{\parindent}{\z@}{-1em}{\subsubparagraphfont}}%
%
\def\@sect#1#2#3#4#5#6[#7]#8{%
\def\tempaa{#8}\def\tempab{}%
\ifx\tempaa\tempab%
\def\@seccntformat##1{\csname the##1\endcsname\hskip4.3pt}%\hskip 1em
\else%
%\def\@seccntformat##1{\csname the##1\endcsname.\hskip1em}%\hskip 1em%
\def\@seccntformat##1{\csname the##1\endcsname\ifappendix\else.\hskip4pt\fi}%\hskip 1em%
\fi%
\ifnum #2>\c@secnumdepth%
     \let\@svsec\@empty\else%
     \refstepcounter{#1}%
     \let\@@protect\protect%
     \def\protect{\noexpand\protect\noexpand}%
     \edef\@svsec{\@seccntformat{#1}}%
     \edef\@ssvsec{\@subseccntformat{#1}}%
     \edef\@prgsec{\@prgformat{#1}}%
     \let\protect\@@protect\fi%
     \@tempskipa #5\relax%
     \ifdim \@tempskipa>\z@%
        \begingroup #6\relax%
     \ifnum#2=3
        \noindent{\hskip #3\relax\@svsec}% {\csname #1numfont\endcsname\@svsec}
                            {\interlinepenalty \@M #8\par}%
       \else
     \ifnum#2=2
        \noindent{\hskip #3\relax\@svsec}% {\csname #1numfont\endcsname\@svsec}
                  {\interlinepenalty \@M #8\par}%
       \else
       \noindent{\hskip #3\relax\@svsec}% {\csname #1numfont\endcsname\@svsec}
                  {\interlinepenalty \@M \uppercase{\textsf{#8}}\par}%
        \fi
      \fi
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname}\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% Use a period (any punc) here for run-in head
                   \@svsec #8.\hskip0pt\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                           \protect\numberline{\csname the#1\endcsname}%
                                     \fi
                       #7}}\fi
     \@xsect{#5}}
%
%%%%%%%%%%%%%%%%%%%%%%%% End Sectioning commands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                LISTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newskip\topsepi \topsepi12\p@ \@plus2\p@% \@minus.5\p@
\newskip\topsepii \topsepii2pt% \@plus1\p@
\newskip\topsepiii \topsepiii2pt% \@plus1\p@
\newskip\itemsepi \itemsepi0pt
\newskip\itemsepii \itemsepii0pt
\newskip\itemsepiii \itemsepiii0pt
\newskip\itemsepiv \itemsepiv0pt
\newskip\itemsepv \itemsepv0pt
\newskip\itemsepvi \itemsepvi0pt
\newskip\itemsepvii \itemsepvii0pt
\newskip\itemsepviii \itemsepviii0pt
\newskip\itemsepix \itemsepix0pt
\newdimen\LabelSep \LabelSep4.7pt
%
\def\@listI{\leftmargin\leftmargini
            \labelwidth\leftmargini
            \advance\labelwidth-\labelsep
            \parsep 0\p@%
            \topsep \topsepi
            \itemsep\itemsepi}%
\let\@listi\@listI
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep\topsepii
              \parsep 0pt
              \itemsep\itemsepii}
\def\@listiii {\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep\topsepiii
              \parsep 0pt
              \itemsep\itemsepiii}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
\def\@listvii {\leftmargin\leftmarginvii
              \labelwidth\leftmarginvii
              \advance\labelwidth-\labelsep}
\def\@listviii {\leftmargin\leftmarginviii
              \labelwidth\leftmarginviii
              \advance\labelwidth-\labelsep}
\def\@listix {\leftmargin\leftmarginvix
              \labelwidth\leftmarginvix
              \advance\labelwidth-\labelsep}
%
\setlength\leftmargini  {2.5em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {1em}
\setlength\leftmarginvi {1em}
\let\leftmarginvii\leftmarginvi
\let\leftmarginviii\leftmarginvi
%%%%%%
%
\setlength  \labelsep  {\LabelSep}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
%
\newcounter{enumv}
\newcounter{enumvi}
\newcounter{enumvii}
\newcounter{enumviii}
\newcounter{enumix}
\renewcommand\theenumi{\arabic{enumi}}
\renewcommand\theenumii{\theenumi.\arabic{enumii}}
\renewcommand\theenumiii{\theenumii.\arabic{enumiii}}
\renewcommand\theenumiv{\theenumiii.\arabic{enumiv}}
\renewcommand\theenumv{\theenumiv.\arabic{enumv}}
\renewcommand\theenumvi{\theenumv.\arabic{enumvi}}
\renewcommand\theenumvii{\theenumvi.\arabic{enumvii}}
\renewcommand\theenumviii{\theenumvii.\arabic{enumviii}}
\renewcommand\theenumix{\theenumviii.\arabic{enumix}}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{\theenumii.}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\newcommand\labelenumv{\theenumv.}
\newcommand\labelenumvi{\theenumvi.}
\newcommand\labelenumvii{\theenumvii.}
\newcommand\labelenumviii{\theenumviii.}
\newcommand\labelenumix{\theenumix.}
%
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumii(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\renewcommand\p@enumv{\p@enumiv\theenumiv}
\renewcommand\p@enumvi{\p@enumv\theenumvi}
\renewcommand\p@enumvii{\p@enumvi\theenumvii}
\renewcommand\p@enumviii{\p@enumvii\theenumviii}
\renewcommand\p@enumix{\p@enumviii\theenumix}


\font\lcirnew = cmsy10 at 10.3pt
\newcommand\bulls{\hskip-2.5pt\raise-.4pt\hbox{\lcirnew\char'017}\hskip1.5pt}


\def\textbullet{\leavevmode\raise2pt\hbox{\hskip2.1pt\bulls}}
\def\textendash{{\bf--}}
\def\textasteriskcentered{\leavevmode\raise-1.5pt\hbox{*}}
\def\textperiodcentered{\leavevmode\raise1.5pt\hbox{\bulls}}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
%
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                \leftmargin10pt        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand\descriptionlabel[1]{\hspace{2pt}%\labelsep
                                \normalfont\bfseries #1}
%
\newenvironment{inlinestripns}[1]{}{}
%
\def\materialsfont{\fontsize{9}{11}\selectfont}
\newenvironment{materialsandmethods}{%
\MaterialsandMethodstrue
    \par%
    \materialsfont}{%
    \par}
%
\def\salutation#1{\gdef\@salutation{#1}}
\salutation{}
\newenvironment{verse}
               {\let\\=\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item[]}
               {\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             ENUMERATE LISTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newskip\listtextleftmargin\listtextleftmargin 0pt%24pt
\newskip\listtextleftmarginii\listtextleftmarginii0pt% 24pt
\newskip\listtextleftmarginiii\listtextleftmarginiii0pt% 24pt
\newskip\listtextrightmargin\listtextrightmargin0pt%.5pc
\newskip\listlabelleftskip \listlabelleftskip0pt%3.3pt
\newskip\listlabelleftskipii \listlabelleftskipii8pt%3.3pt
\newskip\listlabelleftskipiii \listlabelleftskipiii0pt%3.3pt
\newskip\abovelistskipi\abovelistskipi7pt%11pt plus2pt minus1pt
\newskip\belowlistskipi\belowlistskipi7pt%10.5pt plus2pt
\newskip\abovelistskipii\abovelistskipii0pt%6pt plus2pt
\newskip\belowlistskipii\belowlistskipii0pt%6pt plus2pt
\newskip\abovelistskipiii\abovelistskipiii0pt%6pt plus2pt
\newskip\belowlistskipiii\belowlistskipiii0pt%6pt plus2pt
\newskip\abovelistskipiv\abovelistskipiv0pt%6pt plus2pt
\newskip\belowlistskipiv\belowlistskipiv0pt%6pt plus2pt
\newskip\abovelistskipv\abovelistskipv0pt%6pt plus2pt
\newskip\belowlistskipv\belowlistskipv0pt%6pt plus2pt
\newskip\abovelistskipvi\abovelistskipvi0pt%6pt plus2pt
\newskip\belowlistskipvi\belowlistskipvi0pt%6pt plus2pt
\newskip\abovelistskipvii\abovelistskipvii0pt%6pt plus2pt
\newskip\belowlistskipvii\belowlistskipvii0pt%6pt plus2pt
\newskip\abovelistskipviii\abovelistskipviii0pt%6pt plus2pt
\newskip\belowlistskipviii\belowlistskipviii0pt%6pt plus2pt
\newskip\abovelistskipvix\abovelistskipvix0pt%6pt plus2pt
\newskip\belowlistskipix\belowlistskipix0pt%6pt plus2pt
%
%
\newskip\labelsepi \labelsepi4.5pt
\newskip\labelsepii \labelsepii5.2pt
\newskip\labelsepiii \labelsepiii5.2pt%\z@
\newskip\labelsepiv \labelsepiv5.2pt%\z@
\newskip\labelsepv \labelsepv5.2pt%\z@
\newskip\labelsepvi \labelsepvi5.2pt%\z@
\newskip\labelsepvii \labelsepvii5.2pt%%\z@
\newskip\labelsepviii \labelsepviii5.2pt%%\z@
\newskip\labelsepix \labelsepix5.2pt%%\z@
\newdimen\enumdimwd
\newif\iflabelrightalign\labelrightaligntrue
%
\def\listdevicefonti{}
\def\listdevicefontii{}
\def\listdevicefontiii{}
\def\listfont{}
\def\textlistlabel{}
%
\def\enummax#1{%
  \labelsep\csname labelsep\romannumeral\the\@enumdepth\endcsname
  \ifdim\listtextleftmargin>\z@\labelsepi0pt\fi
  \ifdim\listtextleftmarginii>\z@\labelsepii0pt\fi
  \ifdim\listtextleftmarginiii>\z@\labelsepiii0pt\fi
  \setbox\tempbox\hbox{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname#1\hskip\labelsep}%
  \enumdim\wd\tempbox
  \setbox\tempbox\hbox{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname#1}%
  \enumdimwd\wd\tempbox
  \expandafter\global\csname leftmargin\romannumeral\the\@enumdepth\endcsname\enumdim
  \ifdim\listtextleftmargin>\z@
     \leftmargini\listtextleftmargin
     \ifdim\listlabelleftskip>\z@
        \advance\leftmargini-\listlabelleftskip
     \fi
  \fi
  \ifdim\listtextleftmarginii>\z@
     \leftmarginii\listtextleftmarginii
     \ifdim\listlabelleftskipii>\z@
        \advance\leftmarginii-\listlabelleftskipii
     \fi
  \fi
  \ifdim\listtextleftmarginiii>\z@
     \leftmarginiii\listtextleftmarginiii
     \ifdim\listlabelleftskipiii>\z@
        \advance\leftmarginiii-\listlabelleftskipiii
     \fi
  \fi
  \ifdim\listlabelleftskip>\z@
     \advance\leftmargini\listlabelleftskip
  \fi
  \ifdim\listlabelleftskipii>\z@
     \advance\leftmarginii\listlabelleftskipii
  \fi
    \ifdim\listlabelleftskipiii>\z@
     \advance\leftmarginiii\listlabelleftskipiii
  \fi}
%
\enummax{1.}
%
\def\list#1#2{\ifnum \@listdepth >8\relax \@toodeep
     \else \global\advance\@listdepth\@ne \fi
  \rightmargin \z@ \listparindent\z@ \itemindent\z@
  \csname @list\romannumeral\the\@listdepth\endcsname
  \def\@itemlabel{#1}\let\makelabel\@mklab \@nmbrlistfalse #2\relax
  \@trivlist
  \parskip\parsep \parindent\listparindent
  \advance\linewidth -\rightmargin \advance\linewidth -\leftmargin
  \advance\@totalleftmargin \leftmargin
  \parshape \@ne \@totalleftmargin \linewidth
  \ignorespaces}
%
\def\enumerate{\@ifnextchar[{\@enumerate}{\@@enumerate}}%%\@enumerate[\csname label\@enumctr\endcsname]
%
\def\@enumerate[#1]{\par%
\setlength\parskip{0pt}
\def\[{\vspace*{-0pt}\relax\ifmmode\@badmath\else
 \begin{trivlist}%
 \@beginparpenalty\predisplaypenalty
 \@endparpenalty\postdisplaypenalty
 \item[]\leavevmode
 \hbox to\linewidth\bgroup $\m@th\displaystyle
 \hskip\mathindent\bgroup\fi}
%
\def\]{\relax\ifmmode \egroup $\hfil
       \egroup \end{trivlist}\else \@badmath \fi\vspace*{-0pt}}
\global\listingstrue%
      \ifnum \@enumdepth >9 \@toodeep\else%tarun
      \advance\@enumdepth\@ne%
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \setcounter{\@enumctr}{1}\enummax{#1}%
      \list{\csname label\@enumctr\endcsname}{\usecounter{\@enumctr}%
      \topsep\csname abovelistskip\romannumeral\the\@enumdepth\endcsname%
      \itemsep\csname itemsep\romannumeral\the\@enumdepth\endcsname%
      \listfont%\listparindent18.25pt
      \listparindent1em
      \ifnum \@enumdepth=1
         \rightmargin\listtextrightmargin
         \advance\rightmargin\rightskip
         \advance\leftmargin\leftskip
         \tempdimen\leftmargini
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskip{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%\textlistlabel
      \else
      \ifnum \@enumdepth=2
         \tempdimen\leftmarginii
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskipii{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%
      \else
         \ifnum \@enumdepth=3
         \tempdimen\leftmarginiii
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskipiii{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%
     \else
         \def\makelabel##1{\hss\llap{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname##1}}%
      \fi\fi
      \fi}\fi\setlength\parskip{1pt}}
%
\def\endenumerate{\@topsepadd\csname belowlistskip\romannumeral\the\@enumdepth\endcsname\endlist}%
%
\def\@@enumerate{%
\message{^^J^^JERROR:^^J}
\message{****************************************************}
\message{List Environment Must have a Square Bracket Argument}
\message{****************************************************^^J}
\message{^^J}
\message{Example:^^J}
\message{\string\begin{enumerate}[12.]}
\message{^^J----------------------------------------------------^^J}
\message{^^JPress [ENTER] and [e] to correct the error}
\message{^^J-->}
\read-1 to \messageskip}
%
\newenvironment{orderedlist}{%
    \ifnum\@enumdepth=0%
        \mathindent\leftmargini%0pt
        \emathindent\leftmargini%
\else%
    \ifnum\@enumdepth=1%
        \mathindent\leftmarginii%0pt
        \advance\mathindent\labelwidth
        \advance\mathindent\labelsepi
        \emathindent\leftmarginii%
        \advance\emathindent\labelwidth%
        \advance\emathindent\labelsepi%
\else%
    \ifnum\@enumdepth=2%
        \mathindent\leftmarginiii%0pt
        \advance\mathindent\enumdim%
        \advance\mathindent\labelwidth%
        \advance\mathindent-2.5pt
        \emathindent\leftmarginiii%
        \advance\emathindent\labelwidth%
        \advance\emathindent\enumdim%
        \advance\emathindent-2.5pt
\fi\fi\fi%
\overfullrule0pt
\abovelistskipi6pt plus2pt minus.2pt
\belowlistskipi6pt plus2pt minus.2pt
%\listtextrightmargin0pt
\labelsepi5pt
\labelsepii5.3pt
%\listlabelleftskip12pt
  \def\theenumi{\arabic{enumi}}
  \def\theenumii{\theenumi.\arabic{enumii}}%
  \def\theenumiii{\theenumii.\arabic{enumiii}}%
  \def\theenumiv{\theenumiii.\arabic{enumiv}}%
  \def\theenumv{\theenumiv.\arabic{enumv}}%
  \def\theenumvi{\theenumv.\arabic{enumv}}%
  \def\theenumvii{\theenumvi.\arabic{enumv}}%
  \def\theenumviii{\theenumvi.\arabic{enumv}}%
  \renewcommand\labelenumi{\theenumi}
  \renewcommand\labelenumii{\theenumii}
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{wherelist}{%
  \def\listdevicefonti{\upshape}
  \listtextleftmargin 0pt
  \abovelistskipi13\p@ plus2pt
  \belowlistskipi13\p@ plus2pt
  \begin{enumerate}%
}{%
  \end{enumerate}}
  %
\newenvironment{romanlist}{%
  \def\theenumi{\roman{enumi}}\def\theenumii{\roman{enumii}}%
  \def\theenumiii{\roman{enumiii}}\def\theenumiv{\roman{enumiv}}%
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{alphalist}{%
  \def\theenumi{\alph{enumi}}\def\theenumii{\alph{enumii}}%
  \def\theenumiii{\alph{enumiii}}\def\theenumiv{\alph{enumiv}}%
  \renewcommand\labelenumi{(\theenumi)}
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{Romanlist}{%
  \def\theenumi{\Roman{enumi}}\def\theenumii{\Roman{enumii}}%
  \def\theenumiii{\Roman{enumiii}}\def\theenumiv{\Roman{enumiv}}%
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{Alphalist}{%
  \def\theenumi{\Alph{enumi}}\def\theenumii{\Alph{enumii}}%
  \def\theenumiii{\Alph{enumiii}}\def\theenumiv{\Alph{enumiv}}%
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unnumlist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newif\ifboxenv
\global\boxenvfalse
%
\newenvironment{unnumlist}{%
    \list{}{%
%    \listparindent18pt\raggedright
        \leftmargin18pt
        \advance\leftmargin1em
        \itemindent-1em
    \parsep0pt
    \itemsep0pt
    \partopsep0pt}
    \def\makelable##1{##1}%
}{%
    \endlist
}%
%
\def\webadd#1{\par\addvspace{3pt}{\fontsize{9}{12}\selectfont\hskip18pt#1}}
%
\newenvironment{tablist}{%
    \list{}{%
%    \listparindent18pt\raggedright
        \leftmargin18pt
        \advance\leftmargin1em
        \itemindent-1em
    \parsep0pt
    \itemsep0pt
    \partopsep0pt}
    \def\makelable##1{##1}%
}{%
    \endlist
}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Itemize (bullet) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% define \itemizeiTopSep if it is different from \topsepi
%
\newdimen\itemleftmargin \itemleftmargin\z@ % won't be active in enumerate
\newdimen\itemlabelsep \itemlabelsep3.5\p@ % labelsep in itemize for all levels
%
\def\itemize{%
\setlength\parskip{0pt}
  \ifnum \@itemdepth >3 \@toodeep\else\ifnum\@enumdepth>0\@itemdepth\@enumdepth\fi \advance\@itemdepth \@ne
     \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
     \list{\csname\@itemitem\endcsname}%
        {%
        \topsep\csname abovelistskip\romannumeral\the\@itemdepth\endcsname%\fi
        \itemsep\csname itemsep\romannumeral\the\@itemdepth\endcsname
        \labelsep\itemlabelsep
        \listfont
        \setbox\tempbox\hbox{\csname\@itemitem\endcsname}
        \csname leftmargin\romannumeral\the\@itemdepth\endcsname\wd\tempbox
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemlabelsep
        \ifnum\@itemdepth=\@ne
        \ifnum\@enumdepth=0 % if not inside an enumerate
           \rightmargin\listtextrightmargin
           \advance\rightmargin\rightskip
           \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmargini\leftskip
                \else
                \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmarginii\leftskip
        \fi\fi
        \leftmargin\csname leftmargin\romannumeral\the\@itemdepth\endcsname
                  \def\makelabel##1{\hss\llap{##1}}}
   \fi\setlength\parskip{1pt}}%
%
\def\enditemize{%
\@topsepadd\csname belowlistskip\romannumeral\the\@itemdepth\endcsname%\fi
\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Enumroman (i) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{xbllist}{%
   \abovelistskipii9pt plus2pt%
   \belowlistskipii9pt plus2pt%
   \itemlabelsep12.5\p@%
   \itemsepi0pt
   \def\listfont{}
   \renewcommand\labelitemi{\textbullet}%
   \renewcommand\labelitemii{\textbullet}%
   \renewcommand\labelitemiii{\textbullet}%
   \renewcommand\labelitemiv{\textbullet}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{bulletlist}{%
   \abovelistskipii9pt plus2pt%
   \belowlistskipii9pt plus2pt%
   \itemlabelsep12.5\p@%
   \itemsepi0pt
   \renewcommand\labelitemi{\textbullet}%
   \renewcommand\labelitemii{\textbullet}%
   \renewcommand\labelitemiii{\textbullet}%
   \renewcommand\labelitemiv{\textbullet}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{hyphenlist}{%
   \renewcommand\labelitemi{\normalfont\bfseries \textendash}%
   \renewcommand\labelitemii{\normalfont\bfseries \textendash}%
   \renewcommand\labelitemiii{\normalfont\bfseries \textendash}%
   \renewcommand\labelitemiv{\normalfont\bfseries \textendash}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{asterisklist}{%
   \renewcommand\labelitemi{\textasteriskcentered}%
   \renewcommand\labelitemii{\textasteriskcentered}%
   \renewcommand\labelitemiii{\textasteriskcentered}%
   \renewcommand\labelitemiv{\textasteriskcentered}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{caselist}{%
   \def\textlistlabel{Case\ }%
   \begin{enumerate}
}{%
   \end{enumerate}
}
%
\newenvironment{steplist}{%
   \def\textlistlabel{Step\ }%
   \begin{enumerate}
}{%
   \end{enumerate}
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% npara %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\NPARfont{}
\newcounter{ncount}
\newenvironment{npara}{\par%\addvspace{9pt plus2pt}%
   \NPARfont
%   \hsize\textwidth%27pc
   \setcounter{ncount}{0}%
   \def\item{\par\addtocounter{ncount}{1}\arabic{ncount}.\hskip2.6pt\ignorespaces}%
}{%
   \par}%\addvspace{9pt plus2pt}\@endparenv}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unorderedlist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\unorderedlist{\@ifnextchar[{\@unorderedlist}{\@unorderedlist[]}}
\def\@unorderedlist[#1]{\par%
\mathindent11pt
\emathindent11pt
\def\hone{#1}
\ifx\hone\@empty
\leftskip5.6pt
\renewcommand\theenumi{\bulls}
\renewcommand\labelenumi{\theenumi}
\else
\fi
\overfullrule0pt
\setlength{\parindent}{0pt}
\abovelistskipi6pt plus2pt minus.2pt
\belowlistskipi6pt plus2pt minus.2pt
\labelsepi5pt
\labelsepiii5pt
\labelsepiv5pt
\labelsepv5pt
\labelsepvi5pt
\labelsepvii5pt
\labelsepviii5pt
\begin{enumerate}[#1]
\leftskip2pt}
\def\endunorderedlist{\end{enumerate}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\marginfont{\fontsize{9}{10}\selectfont\bfseries\raggedright}
\newcommand{\comment}[1]{%
 \setlength\overfullrule{0pt}%
 \ifodd\c@page
    \marginpar{%
    \parbox{8pc}{\par\vskip-1.5pt\marginfont please check.\par}}%
    \fboxsep.1pt\fbox{#1}%
 \else
    \marginpar{%
    \parbox{8pc}{\par\vskip-1.5pt\marginfont please check.\par}}
    \fboxsep.1pt\fbox{#1}\fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                FLOATS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\setcounter{topnumber}{5}
\def\topfraction{0.99}
\def\textfraction{0.05}
\def\floatpagefraction{0.9}
\setcounter{bottomnumber}{5}
\def\bottomfraction{0.99}
\setcounter{totalnumber}{10}
\def\dbltopfraction{0.99}
\def\dblfloatpagefraction{0.9}
\setcounter{dbltopnumber}{5}

\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{15\p@}
\setlength\belowcaptionskip{3.5\p@}
\newskip\belowfigureskip
\belowfigureskip=4.5pt plus3pt minus1pt
\newskip\abovefigureskip
\abovefigureskip=4.5pt plus3pt minus1pt
%
\def\FigName{figure}
%
\long\def\@makecaption#1#2{\@bsphack%
    \ifx\FigName\@captype%
      \vskip\abovecaptionskip%
    \@makefigurecaption{#1}{#2}%
%\vskip\belowfigureskip
    \else%
      \@maketablecaption{#1}{#2}%
      \vskip\belowcaptionskip%
    \fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\figno#1{\gdef\@figno{#1.\hskip2pt}}
\figno{}
\newcounter{figure}[chapter]
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{{\sc{Figure}}~\thefigure:\ }%
\newenvironment{figure}
               {\@bsphack\@float{figure}}
               {\end@float\@esphack}
\newenvironment{figure*}
               {%\global\ThreadingNotNeededtrue
                \@bsphack\@dblfloat{figure}}
               {\end@dblfloat\@esphack}%\global\ThreadingNotNeededfalse}
%
\long\def\sidepagecaption#1{\insert\footins{%
\refstepcounter{figure}
\vskip4pt
\ifodd\c@page
\noindent{\figcaptionfont{}\hskip3.4pt#1}%
\else
\hfill\noindent{\figcaptionfont{}\hskip3.4pt#1}%
\fi}}
%
%
\def\caption{\refstepcounter\@captype \@dblarg{\@caption\@captype}}
\long\def\@caption#1[#2]#3{\par\begingroup\@bsphack%
%    \@parboxrestore%
    \normalsize%
    \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}\par%
  \endgroup\@esphack}
%
\newdimen\figcapdim
\newbox\onebox

\def\@makefigurecaption#1#2{\figcaptionfont\tbzdbondtrue%
    \setbox\tempbox\hbox{\figcaptionfont{\figcaptionnumfont#1}#2\par}%
\figcapdim\wd\tempbox%
     \ifdim\figcapdim >\hsize%
\setbox\onebox\vbox{\leftskip0pt plus.5fill\rightskip0pt plus-.5fill{\figcaptionfont{\figcaptionnumfont#1} #2\vphantom{y}\parfillskip0pt plus 1fill\par}}%
{\centerline{\box\onebox}}
\else
        \vskip-2.8pt\centerline{\figcaptionfont{\figcaptionnumfont#1} #2}%
\fi\global\tbzdbondfalse}
%
\newbox\figtempbox
%
\def\figurebox#1#2#3{\@ifnextchar[{\@figurebox{#1}{#2}{#3}}{\@figurebox{#1}{#2}{#3}[]}}
\def\TIFFfigurebox#1#2#3{\@ifnextchar[{\@TIFFfigurebox{#1}{#2}{#3}}{\@TIFFfigurebox{#1}{#2}{#3}[]}}
%
\ifIgnoreFigure
    \def\@figurebox#1#2#3[#4]{%
    \MakeRGBcolor[Green]{%
    \hbox to \hsize{\hfill\fboxrule2pt\fbox{\MakeRGBcolor[Red]{\vbox to 4pc{\vfill\hbox to 10pc{\hfill\vbox{\hbox to 10pc{\hfill#4\hfill}\hbox to 10pc{\hfill to be placed here\hfill}\hbox to 10pc{\hfill(Turn off \MakeRGBcolor[Blue]{{\bf ignorefigure}}\hfill}\hbox to 10pc{\hfill option before pagination)\hfill}}\hfill}\vfill}}}%
    \hfill}}}%
\else%
    \def\@figurebox#1#2#3[#4]{%
            \gdef\@thirdarg{#3}%
                \gdef\@frtharg{#4}%
          \ifx\@frtharg\empty%
           \global\figheight=#1%
             \global\figwidth=#2%
          \else%
                     \ifx\@thirdarg\empty\else\epsfxsize#3\fi%
           \setbox\figtempbox=\hbox{\epsfbox{#4}}%
             \global\figwidth=\wd\figtempbox%
           \global\figheight=\ht\figtempbox%
                \fi%
          \vbox{\centerline{\figbox}}}%
\fi%
%\fi
%
    \def\@TIFFfigurebox#1#2#3[#4]{%
            \gdef\@thirdarg{#3}%
                \gdef\@frtharg{#4}%
          \ifx\@frtharg\empty%
           \global\figheight=#1%
             \global\figwidth=#2%
          \else%
                \fi%
          \vbox{\centerline{\inserttiff{#1}{#2}{#4}}}
}%
%
\def\figbox{%
     \ifx\@frtharg\empty%
       \noindent\vbox{\hsize\figwidth%
          \hrule\hbox to\figwidth{\vrule\hfill\vbox to\figheight{\hsize\figwidth\vfill}\vrule}\hrule}%
     \else%
        \vbox{\vskip.8pt\hsize\figwidth%
         \hbox to\figwidth{\vbox to \figheight{\hsize\figwidth\box\figtempbox}}}%
     \fi%
}

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\TBLNO#1{\gdef\@tblno{\fontsize{9}{11}\selectfont{{#1}:\hskip5pt}}}\TBLNO{}%%
%\tblno{}

\font\opf=cmsy6
\def\Tdag{\hbox{\opf\char'171}}


%\def\tref#1{\let\dag\Tdag\ifmmode\hyperlink{#1}{\hbox{$^{\rm #1}$}}\else\hyperlink{#1}{$^{\rm #1}$}\fi}
\def\tref#1{\let\dag\Tdag\ifmmode\hyperlink{\thetabfoot#1}{\hbox{$^{\rm #1}$}}\else\hyperlink{\thetabfoot#1}{$^{\rm #1}$}\fi}
\def\legend#1{\par\vskip.6pt{\fontsize{8}{10}\selectfont\advance\tableleftskip5pt\advance\tablerightskip5pt\leftskip\tableleftskip\rightskip\tablerightskip #1\par}}

\def\textlegend#1{\vskip-4pt\par\noindent{\tabnotefont #1\par}}

\newcounter{tabfoot}%
\newcounter{table}[chapter]
\renewcommand\thetable{\@arabic\c@table}%\thechapter.
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}

\def\fnum@table{{\tablename}~\thetable:\hskip2.6pt}
\newenvironment{table}
               {\addtocounter{tabfoot}{1}\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\addtocounter{tabfoot}{1}%\global\ThreadingNotNeededtrue
               \@dblfloat{table}}
               {\end@dblfloat}%\global\ThreadingNotNeededfalse}
\def\@maketablecaption#1#2{%
    {\tablecaptionfont\noindent{\tablecaptionnumfont #1}%#1
    #2\vphantom{y}\par}}%
%
\newbox\tabbodybox
\newdimen\tabbodyht
\newdimen\tabbodywd
\newskip\askip
\askip15pt

\def\tbl#1#2{\tbzdbondtrue%
\setlength\parskip{0pt}
\def\@xarraycr{\@ifnextchar[\@argarraycr{\ifnum0=`{\fi}${}\cr\noalign{\vskip0pt}}}%
\def\figurebox##1##2##3[##4]{\epsfbox{##4}}
  \global\setbox\tempbox\hbox{\tablefont#2}%
  \tabbodyht\ht\tempbox
  \tabbodywd\wd\tempbox
  \multiply\tabbodyht2
  \advance\tabbodyht-9.5pt
  \tabledim\hsize
  \advance\tabledim-\wd\tempbox
  \ifdim\tabledim>0pt
  \divide\tabledim2
  \else
  \global\tabledim6pt
  \fi
  \global\tableleftskip\tabledim
  \global\tablerightskip\tabledim
  \def\tblcap{#1}
  \caption{\tblcap}%
\ifx\tblcap\empty
\vskip-10pt
\fi
\vskip-10pt
  {\centerline{\tablefont\box\tempbox\vspace*{3pt}}}%
  \global\tbzdbondfalse\setlength\parskip{1pt}}
%
\def\yellowrule{\@ifnextchar[{\@yellowrule}{\@yellowrule[20.6pt]}}
\def\@yellowrule[#1]{\par\hbox to 0pt{\smash{\vrule width\textwidth height0pt depth#1}}}
%
%
\newdimen\tabledim
\newbox\tablecapbox
\newif\iffootnote


\newdimen\tabledim
\newbox\tablecapbox
\newif\iffootnote


\def\TCH#1{\TCHfont#1}%
%
\def\x{@{\extracolsep{\fill}}}
\def\toprule{\\[6.8pt]\Hline\\[-7.4pt]}
\def\colrule{\\[-10.7pt]\Hline\\[-7.3pt]}
\def\botrule{\\[-10.5pt]\Hline\\[-10.7pt]}
\def\crule#1{\\[-9pt]\cline{#1}\\[-6.5pt]}

\def\Hline{%
  \noalign{\ifnum0=`}\fi\hrule \@height .5pt \futurelet%\@height \arrayrulewidth
   \@tempa\@xhline}
%
\def\tabnote#1{{\global\tbzdbondtrue\tabnotefont#1\par}\global\tbzdbondfalse}
\def\tabentry#1{{\tabnotefont\leavevmode\hbox to 12pt{\hfill$^{\rmbox{\rm #1}}$\hskip3pt}}}

\def\tblfno#1{\par\hskip5pt\let\dag\Tdag\fontsize{8}{10}\selectfont\advance\tableleftskip0pt\advance\tablerightskip5pt\leftskip\tableleftskip\rightskip\tablerightskip\hypertarget{\thetabfoot#1}{}\noindent\leavevmode{$^{\rm \ensuremath{#1}}$\hskip0pt}} %manju
\renewenvironment{tabnote}{%
\vskip4.8pt
  \par\parindent9pt\vskip-\lastskip\vskip1.5pt% \ifdim\lastskip=0pt\vskip5pt\else\vskip-\lastskip\fi
  \tablefont
}{%
  \par}
%
\def\tabentry#1{{\tabnotefont$^{\rm #1}$\hskip3pt}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Text Box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcounter{textboxfloat}[chapter]
\def\ftype@textboxfloat{1}
\newif\iftextboxisfloat
\global\textboxisfloatfalse

%%%%%%%%%%%%
\def\fps@textbox{tbp}
\def\ftype@textbox{1}
\def\ext@textbox{lof}
%
\newcounter{boxcount}
\addtocounter{boxcount}{0}
%
\makeatletter
\newif\ifListsp
\global\Listspfalse
% for defining the quote in textbox
\def\textquote{\list{}{\itemindent\z@\topsep6pt
  \leftmargin 22pt \rightmargin10pt}\item[]}
\let\endtextquote=\endlist
%
% for defining the list in textbox
\def\textitemize{\addvspace{6pt}%
  \ifnum \@itemdepth >3
    \@toodeep
  \else
    \advance\@itemdepth \@ne
    \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
    \leftmargini19.75pt%
    \list{\csname\@itemitem\endcsname}%
       {\let\makelabel\right@label}
  \fi\rightskip10pt\global\Listspfalse}%
%
%\let\endtextitemize=\endlist
\def\endtextitemize{\endlist}
\def\textitemize{%
  \ifnum \@itemdepth >3 \@toodeep\else\ifnum\@enumdepth>0\@itemdepth\@enumdepth\fi \advance\@itemdepth \@ne
     \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
     \list{\csname\@itemitem\endcsname}%
        {%
        \topsep\csname abovelistskip\romannumeral\the\@itemdepth\endcsname%\fi
        \itemsep\csname itemsep\romannumeral\the\@itemdepth\endcsname
        \labelsep\itemlabelsep
        \listfont
        \setbox\tempbox\hbox{\csname\@itemitem\endcsname}
        \csname leftmargin\romannumeral\the\@itemdepth\endcsname\wd\tempbox
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemlabelsep
        \ifnum\@itemdepth=\@ne\ifnum\@enumdepth=0 % if not inside an enumerate
           \rightmargin\listtextrightmargin
           \advance\rightmargin\rightskip
           \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmargini\leftskip
        \else
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmarginii\leftskip
        \fi\fi
        \leftmargin\csname leftmargin\romannumeral\the\@itemdepth\endcsname
        \def\makelabel##1{\hss\llap{##1}}}
   \fi}%
%
\def\endtextitemize{%
\@topsepadd\csname belowlistskip\romannumeral\the\@itemdepth\endcsname%\fi
\endlist}

%
\long\def\@makeboxcaption#1#2{%
\parindent0pt
\def\fnum@figure{}
 \setbox\@tempboxa\hbox{{\bfseries Box~\theboxcount.}~#1 #2}%          % previously was commented TPI
    \leftskip10pt\rightskip10pt\unhbox\@tempboxa%   set as justified paragraph
\addtocounter{figure}{-1}}

\global\newdimen\figtopskip
\global\figtopskip0pt
\global\newdimen\figbottomskip
\global\figbottomskip0pt

\newenvironment{textbox}[1][\relax]{%\par%
\leftskip10pt\rightskip10pt%
\let\@makefigurecaption\@makeboxcaption
\ifx#1\relax\begin{figure}\else\begin{figure}[#1]\fi
\global\vspace*{\figtopskip}
\vskip3.5pt
\let\quote\textquote
\let\endquote\endtextquote
\let\itemize\textitemize
\let\enditemize\endtextitemize
%\setlength\parindent{1em}%
\refstepcounter{boxcount}%
\newdimen\tempdimen
\newbox\tempabox
\setbox\tempbox\vbox\bgroup%
\let\it\itshape
\let\bf\bfseries
\selectfont\noindent\vrule height0.5pt width\columnwidth%\par%%
\vskip4.5pt%%
\normalsize%
\leftskip10pt\rightskip10pt%
\ignorespaces}{%
\vskip9pt%
\egroup%
\tempdimen\ht\tempbox%
\advance\tempdimen1.3pt%
\box\tempbox%
\vskip-\tempdimen%
\noindent\vrule height\tempdimen width0.5pt\par\vskip-\tempdimen\vskip-.65pt%
\noindent\vrule height0.5pt width\columnwidth%\hskip-0pt%
\noindent\vrule height\tempdimen width0.5pt\par%
\global\vspace*{\figbottomskip}\end{figure}%
\global\figtopskip0pt\global\figbottomskip0pt}%

\newenvironment{textbox*}[1][\relax]{%\par%
\leftskip10pt\rightskip10pt%
\let\@makefigurecaption\@makeboxcaption
\ifx#1\relax\begin{figure*}\else\begin{figure*}[#1]\fi
\vskip3.5pt
\let\quote\textquote
\let\endquote\endtextquote
\let\itemize\textitemize
\let\enditemize\endtextitemize
\setlength\parindent{10pt}%
\refstepcounter{boxcount}%
\newdimen\tempdimen
\newbox\tempabox
\setbox\tempbox\vbox\bgroup%
\let\it\itshape
\let\bf\bfseries
\selectfont\vrule height0.5pt width\textwidth%\par%%
\vskip4.5pt%%
\normalsize%
\leftskip10pt\rightskip10pt%
\ignorespaces}{%
\vskip7pt%
\egroup%
\tempdimen\ht\tempbox%
\advance\tempdimen1.3pt%
\box\tempbox%
\vskip-\tempdimen%
\noindent\vrule height\tempdimen width0.5pt\par\vskip-\tempdimen\vskip-.65pt%
\noindent\vrule height0.5pt width\textwidth%\hskip-0pt%
\noindent\vrule height\tempdimen width0.5pt\par%
\end{figure*}}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  MATH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\arraycolsep{1.5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
%
\setlength\columnsep{2pc}
\setlength\columnseprule{0\p@}
%
\@addtoreset{equation}{chapter}
\renewcommand\theequation{\arabic{equation}}
\def\@eqnnum{{\reset@font\rmfamily(\theequation)}}
%\def\@eqnnum{{\reset@font\rmfamily()}}
%
\def\bstrut{\vrule width0pt depth6pt}
\def\tstrut{\vrule width0pt height9pt}
\jot=6pt
%
\def\ds{\displaystyle}
\newif\ifcenteredequations
%\global\centeredequationsfalse
\global\centeredequationstrue
\newif\iflistings
\global\listingsfalse

%\def\tag#1{\gdef\@eqnnum{\smash{\normalfont \normalcolor#1}}}

\ifcenteredequations
\else
\mathindent2pc%23.4pt%24pt%
\emathindent2pc%0pt%%
\def\[{%
\@beginparpenalty \predisplaypenalty%
\@endparpenalty   \postdisplaypenalty%
\vskip\abovedisplayskip%
\hbox to \linewidth%
\bgroup$\displaystyle\hskip\mathindent}%
%
\def\]{\hfill$\egroup%
%\par\addvspace{6pt plus2pt minus1pt}%
\vskip\belowdisplayskip%
\@endparenv\noindent\ignorespaces}%
%
\def\[{\vspace*{-6pt}
\relax\ifmmode\@badmath\else
 \begin{trivlist}%
 \@beginparpenalty\predisplaypenalty
 \@endparpenalty\postdisplaypenalty
 \item[]\leavevmode
 \hbox to\linewidth\bgroup $\m@th\displaystyle
 \hskip\mathindent\bgroup\fi}
%
\def\]{\relax\ifmmode \egroup $\hfil
       \egroup \end{trivlist}\else \@badmath \fi
       \vspace*{-6pt}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\renewenvironment{equation}{%
\@beginparpenalty \predisplaypenalty%
\@endparpenalty   \postdisplaypenalty%
\refstepcounter{equation}%
\trivlist \item[]%
\vskip\abovedisplayskip
\ifequation
\vskip-12pt
\else
\vskip0pt
\fi
\hbox to \linewidth%
\bgroup$\displaystyle\hskip\mathindent
}{%
\hfill$\@eqnnum\hskip-\emathindent
\egroup%
\vskip\belowdisplayskip\@endparenv}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\let\eqntopsep\abovedisplayskip
%
\def\equation{\@beginparpenalty\predisplaypenalty
  \@endparpenalty\postdisplaypenalty
\refstepcounter{equation}\trivlist
   \@topsep \eqntopsep                      % used by first \item
   \@topsepadd \eqntopsep                   % used by \@endparenv
   \item[]\leavevmode
  \hbox to\linewidth\bgroup$ \displaystyle \hskip\mathindent\bgroup}
%
\def\endequation{\egroup$\hfil \displaywidth\linewidth
  \@eqnnum\egroup \endtrivlist \addvspace{-\parskip}% SP 1994-12-20
}
%
\def\frac#1#2{{\begingroup#1\endgroup\over#2}}

\def\eqnarray{%\par\vskip-\baselineskip
\stepcounter{equation}%
\def\@currentlabel{\p@equation\theequation}%
\global\@eqnswtrue\m@th%
\global\@eqcnt\z@\tabskip\mathindent%
\let\\\@eqncr%
$$\everycr{}\halign to\displaywidth\bgroup\hskip\@centering%
  $\displaystyle\tabskip\z@skip{##}$\@eqnsel&\global\@eqcnt\@ne%
  \hskip \tw@\arraycolsep%
\hfil${##}$\hfil%
  &\global\@eqcnt\tw@ \hskip \tw@\arraycolsep%
  $\displaystyle{##}$\hfil%
   \tabskip\@centering&\global\@eqcnt\thr@@%
       \hbox to\z@\bgroup\hss##\egroup\tabskip\z@skip\cr}%
%
\def\endeqnarray{%
      \@@eqncr%
      \egroup%
      \global\advance\c@equation\m@ne%
   $$\global\@ignoretrue}

\fi

\newenvironment{fleqnarray}{%
\arraycolsep0pt
\let\rightfill\hfil\begin{eqnarray}%
}{%
   \end{eqnarray}\@ignoretrue}

\def\@arraycr{%
  ${\ifnum0=`}\fi\@ifstar\@xarraycr\@xarraycr}
\def\@xarraycr{\@ifnextchar[\@argarraycr{\ifnum0=`{\fi}${}\cr\noalign{\vskip3pt}}}
\def\@argarraycr[#1]{%
  \ifnum0=`{\fi}${}\ifdim #1>\z@ \@xargarraycr{#1}\else
   \@yargarraycr{#1}\fi}
\def\@xargarraycr#1{\@tempdima #1\advance\@tempdima \dp \@arstrutbox
   \vrule \@height\z@ \@depth\@tempdima \@width\z@ \cr}
\def\@yargarraycr#1{\cr\noalign{\vskip #1}}

\font\sl=linew10 at 7pt
\def\spslas{\hbox{\sl\char'012}}
\def\casefrac#1#2{\raisebox{2.7pt}{{\fontsize{6}{6}\selectfont#1}}\kern-1.5pt\spslas\kern-1.6pt\raisebox{0pt}{{\fontsize{6}{6}\selectfont#2}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Math %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                FOOTNOTE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%% Updated %%%%%%%%%%%%%%%%%%%

\newif\iffootrule
\global\footruletrue

\newif\ifnorule
\global\norulefalse

\iffootrule
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule height .5pt width\ifnorule0pc\else 6pc\fi
  \kern 6\p@}
\else
\def\footnoterule{}
\fi%

\@addtoreset{footnote}{chapter}
\renewcommand\thefootnote{\fnsymbol{footnote}}
%
\def\@fnsymbol#1{%
\ifcase#1\or *\or \dagger\or \ddagger\or
   \SS\or \|\or \#\or **\or \dagger\dagger\or \ddagger\ddagger
   \or \SS\SS \or \|\hskip-1pt\|\or \#\#\or ***\or
   \dagger\dagger\dagger\or \ddagger\ddagger\ddagger\else\@ctrerr\fi\relax}


\def\@fnsymbol#1{\ensuremath{\ifcase#1\or \sstr\or \fdagger\or \fddagger\or
   \fSS  \or \| \or \fpr \or \# \or
 \sstr\sstr\or \fdagger\fdagger   \or \fddagger\fddagger \or \fSS\fSS \or \|\| \or \fpr\fpr\or
 \sstr\sstr\sstr\or \fdagger\fdagger\fdagger   \or \fddagger\fddagger\fddagger \or \fSS\fSS\fSS \or \|\|\| \or \fpr\fpr\fpr\or
\sstr\sstr\sstr\sstr\or \fdagger\fdagger\fdagger\fdagger   \or \fddagger\fddagger\fddagger\fddagger \or \fSS\fSS\fSS\fSS \or \|\|\|\| \or \fpr\fpr\fpr\fpr\or
 \else\@ctrerr\fi}}
%
\newdimen\@footmax
\def\footmax#1{%
  \setbox\tempbox\hbox{$^{#1}$\hskip8.6pt}%
  \global\@footmax\wd\tempbox}
%
\newcount\fnsm
\fnsm0

\long\def\footnote{\@ifnextchar[\@xfootnote{\stepcounter\@mpfn
     \protected@xdef\@thefnmark{\thempfn}%
     \advance\fnsm1%
%     \hyperlink{FNT\the\fnsm}
     {\@footnotemark}\@footnotetext}}

\def\@makefnmark{{\@textsuperscript{\normalfont\@thefnmark}}}

\long\def\@makefntext#1{\par\let\[\(
\let\]\)%
\overfullrule0pt%
    \setbox\tempbox\hbox{#1}%
    \noindent\hskip0pt\hbox{{%\hypertarget{FNT\the\fnsm}{}%%
   \@makefnmark}\hskip1pt}{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Footnote %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\usepackage{multicol}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                PAGE STYLES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\oddfolio{{\rhfont\thepage}}
\def\evenfolio{{\rhfont\thepage}}
%
\def\rightmark{\expandafter\@leftmark\botmark{}{}}
\def\leftmark{\expandafter\@rightmark\firstmark{}{}}
%
\def\ps@headings{%
   \def\@oddfoot{\lower13.7pt\hbox{\@DOInumber}\hfill\lower13.7pt\hbox{\piifont\@jname\,/ \@Fjname}}%
   \def\@evenfoot{\lower13.7pt\hbox{\piifont\@jname\,/ \@Fjname}\hfill\lower13.7pt\hbox{\@DOInumber}}%
   \def\@evenhead{\setlength\overfullrule{0pt}%
   \overfullrule0pt%
   \hbox to \textwidth{\hbox to 0pt{\ifrevise\else\rhfont\evenfolio\fi}\hfill{\rhfont\@rhauthor}\hfill\hbox to 0pt{\llap{\ifearlyview\else\rhfont Vol. \@jvol,\ No. \@issue\fi}}}}%
   \def\@oddhead{\setlength\overfullrule{0pt}%
   \overfullrule0pt%
   \hbox to \textwidth{%
   \hbox to \textwidth{\hbox to 0pt{\ifearlyview\else\rhfont\@jyear\fi}{\hfill{\rhfont\@chaptermark}\hfill\hbox to 0pt{\llap{\ifrevise\else\rhfont\evenfolio\fi}}}}}}%
   \let\@mkboth\markboth
   \def\chaptermark##1{\markboth{##1}{##1}}}
%
\def\ps@template{%
    \def\@evenfoot{\hfill{\dropfoliofont\thepage}}%
    \def\@oddfoot{{\dropfoliofont\thepage}\hfill}%
        \def\@evenhead{{\rhfont\leftmark}\hfil\evenfolio}%
        \def\@oddhead{\oddfolio\hfil{\rhfont\rightmark}}%
    \let\@mkboth\markboth
    \def\chaptermark##1{\markboth{##1}{##1}}%
    \def\sectionmark##1{\markright{##1}}}
    %
\ifprinter
\def\ps@empty{%
      \def\@oddfoot{}%
      \let\@evenfoot\@oddfoot
      \def\@evenhead{}%
      \def\@oddhead{}%
      \let\@mkboth\markboth
      \let\chaptermark\@gobble
      \let\sectionmark\@gobble}
\else
\def\ps@empty{%
      \def\@oddfoot{\hfill\lower\blindfoliodrop\hbox{\thepage}\hfill}
      \let\@evenfoot\@oddfoot
      \def\@evenhead{}%
      \def\@oddhead{}%
      \let\@mkboth\markboth
      \let\chaptermark\@gobble
      \let\sectionmark\@gobble%
}
\fi
%
\def\ps@plain{\let\@mkboth\@gobbletwo%
       \def\@evenfoot{\hfill{\dropfoliofont\thepage}}%
       \def\@oddfoot{{\dropfoliofont\thepage }\hfill}%
       \def\@evenhead{}
       \def\@oddhead{}
     }
%
\def\botfirstpage#1{\gdef\@botfirstpage{#1}}
\botfirstpage{}
\def\botlastpage#1{\gdef\@botlastpage{#1}}
\botlastpage{}
%
\def\@botpagerange{%
\ifx\@botfirstpage\@botlastpage
\@botfirstpage
\else
\@botfirstpage--\@botlastpage\fi}
\def\@botoffprint{\@jid\ \@aid\ \@botpagerange}

\def\ps@copyright{\let\@mkboth\@gobbletwo
       \def\@oddfoot{\ifissue%
        \lower14.6pt\vbox{\hbox to \textwidth{\ifx\@checkcopyrightline\@empty\ifx\@checkFrenchcopyrightline\@empty\else\copyright~\@jyear~\@Frenchcopyrightline\fi\else\ifx\@checkFrenchcopyrightline\@empty\copyright~\@jyear~\@copyrightline\else\copyright~\@jyear~\@copyrightline\,/\,\@Frenchcopyrightline\fi\fi\hfill}}
        \else
        \lower26pt\vbox{\hbox to \textwidth{\ifx\@checkcopyrightline\@empty\ifx\@checkFrenchcopyrightline\@empty\else\copyright~\@jyear~\@Frenchcopyrightline\fi\else\ifx\@checkFrenchcopyrightline\@empty\copyright~\@jyear~\@copyrightline\else\copyright~\ifx\@cjyear\@empty\@jyear\else\@cjyear\fi~\@copyrightline\,/\,\@Frenchcopyrightline\fi\fi\hfill}\par\hbox to \textwidth{\hfill\smash{\lower6pt\hbox{\ifx\@Arno\@empty\else\ifearlyview\else\fbox{\@Arno}\fi\fi}}\hfill}}
        \fi}
       \def\@evenfoot{\lower13.7pt\hbox{\ifx\@checkcopyrightline\@empty\ifx\@checkFrenchcopyrightline\@empty\else\copyright~\@jyear~\@Frenchcopyrightline\fi\else\ifx\@checkFrenchcopyrightline\@empty\copyright~\@jyear~\@copyrightline\else\copyright~\ifx\@cjyear\@empty\@jyear\else\@cjyear\fi~\@copyrightline\,/\,\@Frenchcopyrightline\fi\fi}}%
       \def\@evenhead{\smash{\lower-18.5pt\hbox to 0pt{\oddfolio}}\hfill\hbox{\@catchline}}%
       \def\@oddhead{\hbox to 0pt{\@catchline}\hfill\smash{\ifearlyview\lower-8pt\else\lower-18.5pt\fi\hbox{\oddfolio}}}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                BACK MATTER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Acknowledgement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{ack}[1]{\par\vspace*{-5.4pt}%
\renewcommand\section{\@startsection {section}{1}{\z@}{-17pt plus-1pt minus-1pt}{1.7pt}{\sectionfont}}%
\section*{#1}}{%
\par}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newif\ifappendix
\global\appendixfalse
%

\newcounter{appcount}
\setcounter{appcount}{1}%


\newenvironment{appendix}{\par\appendixtrue%
\section*{APPENDIX\vskip-1.5pt}
\def\SecNo##1{~##1\hskip5pt\ignorespaces}
\def\Secno##1{##1.\hskip5pt\ignorespaces}
\ifnum\theappcount>0
\renewcommand\thesection{APPENDIX}
\else
\renewcommand\thesection{APPENDIX \Alph{section}}
\fi
\renewcommand\thesubsection{}
\renewcommand\theequation{\arabic{equation}}
\setcounter{section}{0}
\setcounter{equation}{0}}{%
\par}



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bibliography Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newdimen\bibitemsep
\newdimen\bibparsep
\newdimen\biblabelsep
\newdimen\bibleftmargini

%%%%%%%% Unnumberd References %%%%%%%%%%%
%
%\setlength\bibleftmargin  {0pt}
\setlength\bibitemsep  {0pt}
\setlength\bibparsep  {10pt}

%%%%%%%% Numbered Within Squre Bracket %%%%%%%%%%%
%
\newif\ifnumsqure
\global\numsquretrue
\setlength\biblabelsep  {3.3pt}
\setlength\bibleftmargini  {3.6pt}
%
%%%%%%%% Numbered With Dot %%%%%%%%%%%%%%%

\def\theBibnum{\arabic{Bibnum}}
\newenvironment{furtherreading}[1]{\par%
    \def\emptyhone{[#1]}%
    \global\setbox\bibbox\hbox{#1.}\global\tempdimen\wd\bibbox%
\ifdim\tempdimen<13pt
\advance\tempdimen2pt
\advance\tempdimen-.8pt
\else
\advance\tempdimen1pt
\bibitemindent-4.8pt
\fi
    \small%
    \section*{Further Reading}%
    \ifx\emptyhone\@empty%
            \def\bibitem##1{\noindent\ignorespaces\hangindent12pt\chaptitlefalse}%
        \else%
            \def\bibitem##1{\par\refstepcounter{Bibnum}\noindent\@hangfrom{%
            \hskip\bibitemindent\leavevmode\hbox to \tempdimen{\hfill[\theBibnum]}%
                \hskip\bibitemskip}\chaptitlefalse\ignorespaces}%
    \fi}{\par}

\newif\ifchaptitle
\global\chaptitlefalse
\def\chaptitle{\global\chaptitletrue}

\newenvironment{thebibliography}[1]{}{}%
\usepackage[round]{Cjsref}
\def\NAT@sort{2}



\def\biblabel#1{\gdef\@biblabel{#1}}
\biblabel{}
\newbox\bibbox
\newskip\bibitemskip
\newskip\bibitemindent
\bibitemskip3.5pt
\def\refname#1{\gdef\@refname{#1}}
\refname{BIBLIOGRAPHY}

\newif\ifRefstyleca
\renewenvironment{thebibliography}[1]{%
\renewcommand\section{\@startsection {section}{1}{\z@}{-12.7pt plus-1pt minus-1pt}{3.5pt}{\sectionfont}}%
\setlength\parskip{0pt}
\def\snm##1{{##1}}
    \par\bibliofont%
    \def\emptyhone{#1}%
    \global\setbox\bibbox\hbox{#1}\global\tempdimen\wd\bibbox%
\ifdim\tempdimen<13pt
\advance\tempdimen1pt
\advance\tempdimen-.8pt
\else
\advance\tempdimen1pt
\bibitemindent-2pt
\fi
    \section*{\@refname}%
     \def\@tempa{#1}%
     \def\@empty{}
     \ifx\@tempa\@empty
     \list{}%
           {\labelwidth0pt\labelsep0pt
           \leftmargin1.3pc%\bibleftmargin%12
           \itemindent-1.3pc%\bibleftmargin
           \itemsep0pt %plus.2pt%\bibitemsep
           \parsep2pt
           \usecounter{enumiv}%
           \let\p@enumiv\@empty
           \renewcommand\theenumiv{\arabic{enumiv}}}%
           \def\@biblabel##1{}
      \else
          \ifnumsqure
          \setbox\tempbox\hbox{[\@tempa]}
          \tempdimen\wd\tempbox
         \ifRefstyleca
                    \setbox\tempbox\hbox{\@tempa\ }
                    \tempdimen\wd\tempbox
                    \def\@biblabel##1{{\hbox to \tempdimen{\hfill##1.}\hfil}}
          \fi
          \else
          \setbox\tempbox\hbox{\@tempa}
          \tempdimen\wd\tempbox
          \def\@biblabel##1{{\hbox to \tempdimen{\hfill##1.}\hfil}}
          \fi%
           \list{\@biblabel{\arabic{enumiv}}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
           \labelsep-3.45pt
           \leftmargin\labelsep
           \advance\leftmargin\tempdimen
          \ifnameyear
          \labelwidth\z@ \labelsep\z@
          \leftmargin18pt
          \itemindent-18pt
        \else
          \ifnumsqure
               \advance\leftmargin\bibleftmargini
               \labelsep\biblabelsep
         \fi\fi
           \itemsep0pt %plus.1pt%\bibitemsep
           \parsep\bibparsep
           \usecounter{enumiv}%
           \let\p@enumiv\@empty
           \renewcommand\theenumiv{\arabic{enumiv}}}%
     \fi
     \sloppy\clubpenalty4000\widowpenalty10000\clubpenalty=10000%
\hbadness=10000\vbadness=10000\hfuzz=10001pt\vfuzz=10001pt\sfcode`\.=\@m%
}{\def\@noitemerr
     {\@latex@warning{Empty `thebibliography' environment}}%
     \endlist\setlength\parskip{1pt}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\alphabet#1{{\bf#1}\par\vskip.5\baselineskip}
%
\newenvironment{theindex}{%
      \columnseprule \z@
      \columnsep 12\p@
      \indexfont\markboth{\indexname}{\indexname}
      \tempdimen\normalbaselineskip
      \advance\tempdimen-\baselineskip
      \twocolumn[\@makeschapterhead{\indexname}\vskip\tempdimen]%
      \parindent\z@
      \let\item\@idxitem
}{%
      \par}
\newcommand\@idxitem  {\par\hangindent3em}
\newcommand\subitem   {\par\hangindent3em\hspace*{1em}}
\newcommand\subsubitem{\par\hangindent3em\hspace*{2em}}
\newcommand\indexspace{\par \vskip\baselineskip}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg {2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{2}
%
\newcommand\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength\@tempdima{3em}%
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
       \nobreak
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}
    \endgroup
  \fi}
\newcommand\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}
\newcommand\l@section      {\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\l@subsection   {\@dottedtocline{2}{3.8em}{3.2em}}
\newcommand\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\newcommand\l@paragraph    {\@dottedtocline{4}{10em}{5em}}
\newcommand\l@subparagraph {\@dottedtocline{5}{12em}{6em}}
\let\l@table\l@figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End FM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Character count %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@questionmark{\vbox to 0pt{\vskip13pc\hskip-11pc%5pc
                  \hbox to \textwidth{\special{color push}\gray\fontsize{570}{570}\selectfont ?
                  \special{color pop}\hss}}}
\def\questionmark#1{\xdef\@questionmark{#1}}%
%
\def\@TotalChars{}
\def\TotalChars#1{\xdef\@TotalChars{#1}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File Version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifprinter
  \def\draftnote{}
  \def\comment{}
  \def\Qauthor#1{}
   \def\QtoAuthor#1{}
   \def\QtoClient#1{}
   \def\QtoQC#1{}
   \def\QtoGL#1{}
   \def\qtoa#1{}
\else
\iftrim
  \def\draftnote{}
  \def\comment{}
  \def\Qauthor#1{}
   \def\QtoAuthor#1{}
   \def\QtoClient#1{}
   \def\QtoQC#1{}
   \def\QtoGL#1{}
   \def\qtoa#1{}
\else
\def\draftnote{%\vbox to 0pt{\vskip-12pt
}%
\def\Qauthor#1{\marginpar{{\raggedright\footnotesize\bf #1\endgraf}}}
\fi\fi
\def\Qauthor#1{\setlength\marginparsep{20\p@}\marginpar{{\raggedright\footnotesize #1\endgraf}}}
\let\qtoa\Qauthor
\gdef\@nodocument{}%
%
\def\enddocument{%
   \@enddocumenthook
   \@checkend{document}%
   \clearpage
   \begingroup
     \if@filesw
       \immediate\write\@mainaux{\string\questionmark{}}%
       \immediate\closeout\@mainaux
       \let\@setckpt\@gobbletwo
       \let\@newl@bel\@testdef
       \@tempswafalse
       \makeatletter \input\jobname.aux
     \fi
     \@dofilelist
     \ifdim \font@submax >\fontsubfuzz\relax
       \@font@warning{Size substitutions with differences\MessageBreak
                  up to \font@submax\space have occured.\@gobbletwo}%
     \fi
     \@defaultsubs
     \@refundefined
     \if@filesw
       \ifx \@multiplelabels \relax
         \if@tempswa
           \@latex@warning@no@line{Label(s) may have changed.
               Rerun to get cross-references right}%
         \fi
       \else
         \@multiplelabels
       \fi
     \fi
   \endgroup
   \deadcycles\z@\@@end}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% contact sty %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%§¿¿¿ileif@resent{Pass.Pas}
\def\cropmarks{\nointerlineskip\vbox to 0pt{\offinterlineskip\vss
  \hbox to \trimwidth{\vbox to 0pt{\draftnote}\hfill}}}%

\def\croprules{%
        \nointerlineskip\vbox to 0pt{\offinterlineskip
  \vskip4pt\cropmarks
        \vskip\topmargin
  \hbox to \trimwidth{\hskip\@themargin\hbox to \textwidth{\hss%
  \vrule height0pt depth\typeheight width.1pt%
  \vrule height.1pt depth0pt width\textwidth
  \vrule height0pt depth\typeheight width.1pt\hss}}%
  \hbox to \trimwidth{\hskip\@themargin\hbox to \textwidth{\hss%
                                                                                \vrule height.1pt depth0pt width\textwidth\hss}\hss}
}}
%
\gdef\@placewidthart{}
%
\figwidth0pt
\figheight0pt
%
\def\authorname{}
\def\jobno{}
\def\chapno{}
%
\def\changeh#1.#2pt{}
%
\def\ArtPiece#1#2#3{%
\fboxrule.1pt
\fboxsep0pt
 \oldfigwidth\figwidth
  \gdef\figtemp{#1}%
 \gdef\figtempscale{#3}%\ifx\figtempscale\empty\figtempcount100\else\figtempcount\figtempscale\fi\multiply\figtempcount10
 \setbox\tempbox=\vbox{\epsfbox{\ArtDir#1}}%
    \figheight\ht\tempbox
 \setbox\tempbox=\hbox{\epsfbox{\ArtDir#1}}%
 \advance\figwidth\wd\tempbox
 \presentfigwidth\wd\tempbox
 \remainingtextwidth\textwidth
 \advance\remainingtextwidth-\oldfigwidth
 \advance\remainingtextwidth-\presentfigwidth
 \tempdim\oldfigwidth\advance\tempdim\presentfigwidth
   %\ifdim\tempdim>\textwidth\par\vskip6pt \figwidth0pt
 \ifdim\remainingtextwidth<0pt\par\vskip6pt %\figwidth0pt
 \else
 \fi
 \setbox\tempbox=\hbox{\fbox{\epsfbox{\ArtDir#1}\hss}}%
\xdef\@placewidthart{\@placewidthart\unhbox\tempbox}%
  \leavevmode\vbox{\hsize\presentfigwidth\parindent0pt\leftskip0pt plus1fill\rightskip\leftskip{}
      \@placewidthart\par
      \vskip3pt
%      Fig. #2\par
      Figure: \figtemp\par
      Reduction: \figtempscale\%\par
                        Height:\the\figheight\par
                        Width: \the\presentfigwidth\par
     \jobno\quad \chapno \quad\authorname}\unskip\hskip10pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End File Version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                TRIMMARKS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{linecounter}
\newcommand\thelineno{\arabic{linecounter}}
\setcounter{linecounter}{10000}
%
\def\trimmarks{%
  \vbox to 0pt{%
     \vskip-25pt\parindent0pt
     \draftnote\offinterlineskip}%
  \ifdraft
    {\vbox to 0pt{\hsize\trimwidth\offinterlineskip
    \parindent0pt\leftskip0pt\rightskip0pt
    \fboxsep0pt\fboxrule\draftrule
    \fbox{\vbox{\vskip\trimheight\hskip\trimwidth}}}}%
    {\vbox to 0pt{\hsize\trimwidth\offinterlineskip
    \parindent0pt\leftskip0pt\rightskip0pt
    \vbox to \trimheight{\offinterlineskip\parindent0pt
    \hbox to \trimwidth{\vbox to 2pc{\vskip-3.5pc\hbox{\vrule height2pc width\trimrule}}\raisebox{2pc}{\hbox{\hskip-3.5pc\vrule width2pc height\trimrule}}\hfill
      \raisebox{2pc}{\vrule width2pc height\trimrule\hskip-3.75pc}
      \vbox to 2pc{\vskip-3.5pc\hbox{\vrule height2pc width\trimrule}}%
    }\vfill
    \hbox to \trimwidth{\hbox{\hskip-3.5pc\vrule height\trimrule width2pc}\vbox to 3pc{\vspace*{4.5pc}\hbox{\hskip1.5pc\vrule width\trimrule height2pc}}\hfill
    \vbox to 3.5pc{\vskip5pc\hbox{\vrule height2pc width\trimrule}}\rlap{\hskip1.5pc\vrule width2pc height\trimrule}}}}}%
  \else
    {\vbox to 0pt{\hsize\trimwidth\offinterlineskip
    \parindent0pt\leftskip0pt\rightskip0pt
    \fboxsep0pt\fboxrule\draftrule
    \fbox{\vbox{\vskip\trimheight\hskip\trimwidth}}}}%
    {\vbox to 0pt{\hsize\trimwidth\offinterlineskip
    \parindent0pt\leftskip0pt\rightskip0pt
    \vbox to \trimheight{\offinterlineskip\parindent0pt
    \hbox to \trimwidth{\vbox to 2pc{\vskip-3.5pc\hbox{\vrule height2pc width\trimrule}}\raisebox{2pc}{\hbox{\hskip-3.5pc\vrule width2pc height\trimrule}}\hfill
      \raisebox{2pc}{\vrule width2pc height\trimrule\hskip-3.75pc}
      \vbox to 2pc{\vskip-3.5pc\hbox{\vrule height2pc width\trimrule}}%
    }\vfill
    \hbox to \trimwidth{\hbox{\hskip-3.5pc\vrule height\trimrule width2pc}\vbox to 3pc{\vspace*{4.5pc}\hbox{\hskip1.5pc\vrule width\trimrule height2pc}}\hfill
    \vbox to 3.5pc{\vskip5pc\hbox{\vrule height2pc width\trimrule}}\rlap{\hskip1.5pc\vrule width2pc height\trimrule}}}}}
\fi
\insidedraftrules}
%
\newcount\a
\newcount\b

\newbox\mybox
\newbox\tempbox

\def\GenerateLineNumbers{%\@bsphack%
\a=1
\loop
\ifnum\a<47
\moveleft1.2pc\vbox{
\the\a\vphantom{\big|}}
\advance\a by 1
\repeat
%\@esphack
}

\newdimen\halfcol
\halfcol\textwidth
\divide\halfcol2
\advance\halfcol-.5\columnsep %Tarun

\if@twocolumn
\def\midrule{\vskip\headheight\vskip\headsep\hskip\halfcol\vrule height\normaltextheight width\draftrule\hskip\columnsep\vrule height\normaltextheight width\draftrule}
\else
\def\midrule{}
\fi



\def\insidedraftrules{\setlength\parskip{0pt}\vbox to 0pt{%
  \offinterlineskip\parindent0pt
  \vskip\topmargin
  \tempdimen\normaltextheight
  \advance\tempdimen\headheight
  \advance\tempdimen\headsep
\advance\tempdimen-.1pt
  \moveright\@themargin
  \vbox{\vskip.15pt%
{\vbox to 0pt{%
\vskip\headheight\vskip\headsep% added .15pt for the running head to touch top of draft rule
   \vrule height\draftrule width\textwidth}}
   {\hbox to 0pt{\fboxsep0pt\fboxrule\draftrule
        \fbox{\vbox to \tempdimen
                    {\hsize\textwidth \midrule\hskip\textwidth}}}}
\ifauto
\ifspecial
\smash{\lower30pt\hbox to\textwidth{\small\hfill\@botoffprint}}
\else
\smash{\lower30pt\hbox to\textwidth{\small\hfill\@botoffprint}}
\fi\fi}}\setlength\parskip{1pt}}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Trimmarks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                OUTPUT ROUTINE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def \@shipoutsetup {%
    \@resetactivechars
    \let\-\@dischyph
     \let\'\@acci\let\`\@accii\let\=\@acciii
     \if@specialpage
       \global\@specialpagefalse\@nameuse{ps@\@specialstyle}%
    \fi
     \if@twoside
       \ifodd\count\z@ \let\@thehead\@oddhead \let\@thefoot\@oddfoot
            \let\@themargin\oddsidemargin
       \else \let\@thehead\@evenhead \let\@thefoot\@evenfoot
        \let\@themargin\evensidemargin
       \fi
    \fi
     \reset@font
     \normalsize
     \baselineskip\z@skip \lineskip\z@skip \lineskiplimit\z@
     \let\par\@@par
}
%
\newcount\PageCountForBookmarks
\global\PageCountForBookmarks=1

\def\@outputpage{%
  \let \protect \@unexpandable@protect
  \shipout \vbox{%
    \let \protect \relax
    \aftergroup\let \aftergroup\protect \aftergroup\relax
\stepcounter{linecounter}
    \@shipoutsetup\trimmarks%
    \vskip \topmargin
    \moveright\@themargin \vbox {%
      \setbox\@tempboxa \vbox to\headheight{%
        \vfil
        \color@begingroup
          \normalcolor
          \hbox to\textwidth {%
            \let \label \@gobble
            \let \index \@gobble
            \let \glossary \@gobble %% 21 Jun 91
            \@thehead
            }%
        \color@endgroup
        }%                        %% 22 Feb 87
      \dp\@tempboxa \z@
      \box\@tempboxa
      \vskip \headsep
      \box\@outputbox
      \baselineskip \footskip
      \color@begingroup
        \normalcolor
        \hbox to\textwidth{%
          \let \label \@gobble
          \let \index \@gobble      %% 22 Feb 87
          \let \glossary \@gobble   %% 21 Jun 91
          \@thefoot
          }%
      \color@endgroup
      }%
    }%
  \global \@colht \normaltextheight
  %      \endgroup
  \stepcounter{page}%
\global\PageCountForBookmarks\c@page
\global\advance\PageCountForBookmarks-\@firstpage
\global\advance\PageCountForBookmarks1
  \let\firstmark\botmark
  \ifodd\c@page
     \ifspreadlong\global\spreadlongfalse
     \enlargethispage{\@spreadlong}\fi
  \fi
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Misc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\ignore#1{}
%
\def\@@question#1{\hbox to \hsize{\hfill
                   \rlap{\hskip\marginparsep
                         \vbox to 0pt{\hsize\marginparwidth
                                      \footnotesize
                                      \raggedright#1\endgraf\vss}}}}
\def\@question#1{\ifvmode
                    \@@question{#1}%
                 \else
                     \vadjust{\vbox to 0pt{%
                     \vskip-7.5pt\@@question{#1}\vskip7.5pt}}
                 \fi
}
%
\def\QtoAuthor#1{\@question{{\bf Author:\ #1}}}
\def\QtoClient#1{\@question{{\bf Client/PM:\ #1}}}
\def\QtoQC#1{\@question{{\bf QC:\ #1}}}
\def\QtoGL#1{\@question{{\bf GL:\ #1}}}
\def\qtoR#1{\@question{\hskip12.9pt{#1}}}
\def\QtoA#1{\@question{\vskip-6pt{\bf #1}}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Output Routine  %%%%%%%%%%%%%%%%%%%%%%
%
\newif\ifspreadlong
\def\spreadlong#1{\ifodd\c@page\wlog{Ignoring spreadlong}
                  \else
                     \spreadlongtrue\gdef\@spreadlong{#1}%
                     \enlargethispage{#1}%
                  \fi}
\def\leaflong#1{\enlargethispage{#1}}
%
%%%%%%%%%%%%%%%%%%%%%%%% Proof Readers Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcount\bcount
\newcount\ncount\ncount=1
\newcount\ccount\ccount=1
\newdimen\temptextheight
\newcommand\proofreaderstemplate{\overfullrule0pt%
   \temptextheight\textheight
   \trimrule.25pt
   \markboth{SECTION TITLE}{SECTION TITLE}
   \clearpage
   \thispagestyle{template}   \ptdraftrule.1pt
   \vbox to \z@{\vskip-\headheight\vskip-\headsep\vskip-\topskip%\vskip0.5pt%
   \noindent\hspace*{-\extramargin}\fboxrule\ptdraftrule\fboxsep0pt\fbox{\vbox to 6.5pt{\hbox to\typewidth{}}\hspace*{-2\ptdraftrule}}
   \vskip\headheight\vskip\headsep\vskip-\baselineskip\vskip4.3pt
   \noindent\hspace*{-\extramargin}\fbox{\vbox to \temptextheight{\hbox to\extramargin{}}\hspace*{-2\ptdraftrule}}%
   \fbox{\vbox to \temptextheight{\hbox to\textwidth{}}\hspace*{-2\ptdraftrule}}
   \vskip\footskip
   \vskip-12pt
   \noindent\hspace*{-0pc}\fboxrule\ptdraftrule\fboxsep0pt\fbox{\vbox to 6.5pt{\hbox to\typewidth{}}\hspace*{-2\ptdraftrule}}}%
   \vskip-\baselineskip
      \def\thepage{000}
      \tempdimen\normaltextheight
      \divide\tempdimen\normalbaselineskip
      \bcount\tempdimen
      \advance\bcount1
      \loop
      \ifnum\bcount>0
        \ifnum\ncount=10
          \hspace*{6pc}\the\ccount \vskip0.0000001pt
          \ncount=1
        \else
          \ifnum\ccount>49
             \hspace*{6pc}\the\ccount \vskip0.0000001pt
          \else
             \hspace*{6pc}\the\ncount \vskip0.0000001pt
          \fi
          \advance\ncount1
        \fi
        \advance\bcount-1
        \advance\ccount1
      \repeat
   \ncount=1\ccount=1\bcount=0}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\medmuskip=4mu
\thickmuskip=5mu
\pagestyle{headings}
\pagenumbering{arabic}      % Arabic page numbers
\sloppy
\frenchspacing
%\hoffset=7.6pc
%\voffset=13.5pc
%\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nomenclature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newdimen\ttempht
\newdimen\leftsp
\newdimen\bdimen
\newdimen\botruleht
\newdimen\topruleht
\botruleht.4pt
\topruleht.4pt
\newcount\tempcount
\newcount\bcount
\newbox\ttbox
\newbox\firstbox
\newbox\thirdbox

\def\NomenClatureAboveSpace{14pt plus2pt minus1pt}
\def\NomenBoxSep{10pt}
\def\NomenClatureBelowSpace{14pt plus2pt minus1pt}
\def\NomenClatureTitleAboveSpace{6pt plus2pt minus1pt}
\def\NomenClatureTitleBelowSpace{6pt plus2pt minus1pt}

\def\NomenclatureBaseFont{\leftskip\NomenBoxSep\rightskip\NomenBoxSep plus1fill}
\def\NomenclatureHeadFont{\bfseries}
\def\NomenclatureTitleFont{\itshape}

\def\nomenname#1{\gdef\@nomenname{#1}}
\def\nomenname{Nomenclature}

\def\nomentitle#1{\par\addvspace{\NomenClatureTitleAboveSpace}%
#1\par\addvspace{\NomenClatureTitleBelowSpace}}

\newbox\NomenBox
\newdimen\NomenDim

\def\nomenentry#1#2{\leavevmode\@hangfrom{\hbox to \NomenDim{#1}\quad}#2\hfil\par}

\newenvironment{nomenclature}[1]{\par\addvspace{\NomenClatureAboveSpace}\parindent0pt%
\setbox\NomenBox\hbox{#1}
\NomenDim\wd\NomenBox
   \setbox\tempbox\vbox\bgroup
   \vskip\NomenBoxSep
   \NomenclatureBaseFont
   \noindent{\bf \nomenname}\par\vskip\baselineskip
}{%
   \vskip-\lastskip\vskip\NomenBoxSep
   \egroup
   \ttempht=\ht\tempbox
   \ifdim\pagegoal=\maxdimen
     \leftsp\vsize
   \else
     \leftsp\pagegoal
   \fi
   \advance\leftsp-\pagetotal   \advance\leftsp-18pt
   \wlog{\the\pagegoal}
\ifdim\ttempht < \leftsp
   \botruleht.5pt   \topruleht.5pt
   \makeexamplebox{\ttempht}{\tempbox}%
\else
   \setbox\firstbox=\vsplit\tempbox to \leftsp
   \botruleht0pt   \topruleht.5pt
   \makeexamplebox{\leftsp}{\firstbox}%
   \advance\ttempht by -\textheight   \divide\ttempht by \textheight   \bcount\ttempht
   \loop
   \ifnum\bcount > 0
      \setbox\firstbox=\vsplit\tempbox to \textheight
      \topruleht.0pt
      \makeexamplebox{\textheight}{\firstbox}%
      \advance\bcount by -1
   \repeat
   \ttempht\ht\tempbox
   \botruleht.5pt\topruleht.0pt
   \makeexamplebox{\ttempht}{\tempbox}%
\fi
\addvspace{\NomenClatureBelowSpace}\@endparenv}

\newdimen\nomdim
\renewenvironment{nomenclature}[1]{\par%
\overfullrule0pt%%
\setbox\tempbox\vbox\bgroup%%
\hrule height.5pt width\columnwidth%%
\vskip10pt%%
{\bfseries #1}\par
\leftskip10pt%%
\rightskip10pt%%
}{\par%%
\egroup%%
\nomdim\ht\tempbox
\box\tempbox%%
\par\vskip-\nomdim\vskip-.65pt
\noindent\vrule height\nomdim width.5pt%
\noindent\vrule height.5pt width\columnwidth%\textwidth%\hskip-0pt%
\leavevmode\hspace*{-1pt}\noindent\vrule height\nomdim width.5pt%
\par\addvspace{12pt plus 2pt minus1pt}}

\def\deflist{\@ifnextchar[{\@deflist}{\@deflist[]}}
\newbox\defbox%
\newdimen\defboxdimen%
\def\@deflist[#1]{%
\setbox\defbox\hbox{#1:\quad}%
\defboxdimen\wd\defbox%
\def\tempa{#1}%
\par\addvspace{12pt plus2pt minus2pt}%
\setlength\parindent{0pt}%
\def\deftitle##1{{\noindent\itshape##1}\par}%\fontsize{12}{12}\selectfont
\ifx\tempa\@empty%
\def\defitem##1{\@hangfrom{##1\ }}%%%
\else%
\def\defitem##1{\@hangfrom{\hbox to \defboxdimen{##1\quad}}}%%
\fi%
\def\defterm##1{\gdef\@defterm{##1}}
\def\defdefinition##1{\@hangfrom{\hbox to \defboxdimen{\@defterm\hfill}}{{##1}}\par}
}
\def\enddeflist{\par\addvspace{8pt plus2pt minus2pt}}
%
\def\makeexamplebox#1#2{% #1 ht; #2 box
   \vbox to #1{%
   %\vbox{{\vrule height.5pt width\columnwidth}}
   %\vskip-#1
   \vbox{\hrule height\topruleht\vrule\vbox{\box#2}\kern-.79pt\vrule\hrule height\botruleht}}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           DEFINING NEW ENVIRONMENTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\MakeCommand#1{\csname#1\endcsname}
\def\NEWCounter#1{\newcounter{#1}}
%
\def\NewEnvironment#1{%
\newenvironment{#1}[2]{\par\addvspace{\csname#1@AboveSpace\endcsname}%
\setlength\parindent{\MakeCommand{#1@ParaIndent}}%
\MakeCommand{#1@AboveRule}%
\MakeCommand{#1@StartIndent}%
{\MakeCommand{#1@HeadFont}\MakeCommand{#1@Name}\MakeCommand{#1@HeadToNumSep}##1}%
\def\scndarg{##2}%
\ifx\scndarg\@empty%
{\bfseries\MakeCommand{#1@SecPunctuationStyle}}%
\else%
\MakeCommand{#1@HeadToSectionSep}%%
\MakeCommand{#1@SecOpenStyle}{##2}%%
\MakeCommand{#1@SecCloseStyle}%%
\MakeCommand{#1@SecPunctuationStyle}\fi%}%
%\MakeCommand{#1@Count}%%
%\MakeCommand{#1@Counter}%%
\MakeCommand{#1@HeadParaSep}%%
\MakeCommand{#1@Font}}{%
\MakeCommand{#1@EndLogo}%
\MakeCommand{#1@BelowRule}%
\par\addvspace{\csname#1@BelowSpace\endcsname}%
\@endparenv\reset@font}%
}
%
\def\halmos{$\blacksquare$}   % Set End logo if present

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Proposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{proposition}
\def\proposition@ParaIndent{0pt}
\def\proposition@StartIndent{\leftskip0pt}
\def\proposition@Name{}
\def\proposition@Font{\itshape}
\def\proposition@HeadFont{\bf}
\def\proposition@HeadToNumSep{}
\def\proposition@HeadToSectionSep{\ }
\def\proposition@SecOpenStyle{\bgroup}
\def\proposition@SecCloseStyle{\egroup}
\def\proposition@SecPunctuationStyle{}%{ \thesection .\global\advance\propositioncount by 1
%\the\propositioncount.}
\def\proposition@HeadParaSep{\quad}
\def\proposition@AboveSpace{8pt plus2pt minus2pt}
\def\proposition@BelowSpace{8pt plus2pt minus2pt}


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lemma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{lemma}
\def\lemma@ParaIndent{0pt}
\def\lemma@StartIndent{\leftskip0pt}
\def\lemma@Name{}
\def\lemma@Font{\itshape}
\def\lemma@HeadFont{\bf}
\def\lemma@HeadToNumSep{}
\def\lemma@HeadToSectionSep{\ }
\def\lemma@SecOpenStyle{\bgroup}
\def\lemma@SecCloseStyle{\egroup}
\def\lemma@SecPunctuationStyle{}%{ \thesection .\global\advance\lemmacount by 1
%\the\lemmacount.}
\def\lemma@HeadParaSep{\quad}
\def\lemma@AboveSpace{8pt plus2pt minus2pt}
\def\lemma@BelowSpace{8pt plus2pt minus2pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sublemma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{sublemma}
\def\sublemma@Counter{\csname thesublemma\endcsname}%@Counter{~\thesection.\csname thesublemma\endcsname}
\def\sublemma@ParaIndent{12pt}
\def\sublemma@StartIndent{}
\def\sublemma@Name{Sublemma}
\def\sublemma@Font{\itshape}%\it}
\def\sublemma@HeadFont{\bf}
\def\sublemma@HeadToNumSep{}
\def\sublemma@HeadToSectionSep{\ }
\def\sublemma@SecOpenStyle{\bgroup(}
\def\sublemma@SecCloseStyle{)\egroup}
\def\sublemma@SecPunctuationStyle{.}
\def\sublemma@HeadParaSep{\hskip2pt}
\def\sublemma@AboveSpace{6pt plus2pt minus1pt}
\def\sublemma@BelowSpace{6pt plus2pt minus1pt}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{theorem}
\def\theorem@ParaIndent{0pt}
\def\theorem@StartIndent{\leftskip0pt}
\def\theorem@Name{}
\def\theorem@Font{\itshape}
\def\theorem@HeadFont{\bf}
\def\theorem@HeadToNumSep{}
\def\theorem@HeadToSectionSep{\ }
\def\theorem@SecOpenStyle{\bgroup}
\def\theorem@SecCloseStyle{\egroup}
\def\theorem@SecPunctuationStyle{}%{ \thesection .\global\advance\theoremcount by 1
%\the\theoremcount.}
\def\theorem@HeadParaSep{\quad}
\def\theorem@AboveSpace{8pt plus2pt minus2pt}
\def\theorem@BelowSpace{8pt plus2pt minus2pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{definition}
\def\definition@ParaIndent{0pt}
\def\definition@StartIndent{\leftskip0pt}
\def\definition@Name{}
\def\definition@Font{\itshape}
\def\definition@HeadFont{\bf}
\def\definition@HeadToNumSep{}
\def\definition@HeadToSectionSep{\ }
\def\definition@SecOpenStyle{\bgroup}
\def\definition@SecCloseStyle{\egroup}
\def\definition@SecPunctuationStyle{}%{ \global\advance\definitioncount by 1 \the\definitioncount.}
\def\definition@HeadParaSep{\quad}
\def\definition@AboveSpace{8pt plus2pt minus2pt}
\def\definition@BelowSpace{8pt plus2pt minus2pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Proof %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{proof}
\def\proof@ParaIndent{12pt}
\def\proof@StartIndent{\leftskip0pt}
\def\proof@Name{}
\def\proof@Font{}
\def\proof@HeadFont{\it}
\def\proof@HeadToNumSep{}
\def\proof@HeadToSectionSep{\ }
\def\proof@SecOpenStyle{\bgroup}
\def\proof@SecCloseStyle{\egroup}
\def\proof@SecPunctuationStyle{\rm .}% \thesection .\global\advance\proofcount by 1
%\the\proofcount.
\def\proof@HeadParaSep{\quad}
\def\proof@AboveSpace{8pt plus2pt minus2pt}
\def\proof@BelowSpace{8pt plus2pt minus2pt}
\def\proof@EndLogo{\hspace*{10pt plus1fill}\halmos}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% corollary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{corollary}
\def\corollary@ParaIndent{0pt}
\def\corollary@StartIndent{\leftskip0pt}
\def\corollary@Name{}
\def\corollary@Font{\itshape}
\def\corollary@HeadFont{\bf}
\def\corollary@HeadToNumSep{}
\def\corollary@HeadToSectionSep{\ }
\def\corollary@SecOpenStyle{\bgroup}
\def\corollary@SecCloseStyle{\egroup}
\def\corollary@SecPunctuationStyle{}%{ \thesection .\global\advance\corollarycount by 1
%\the\corollarycount.}
\def\corollary@HeadParaSep{\quad}
\def\corollary@AboveSpace{8pt plus2pt minus2pt}
\def\corollary@BelowSpace{8pt plus2pt minus2pt}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% note %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{Note}
\def\Note@Counter{\csname theNote\endcsname}
\def\Note@ParaIndent{12pt}
\def\Note@StartIndent{}
\def\Note@Name{Note}
\def\Note@Font{\itshape}%\leftskip0pt\rightskip0pt\mathindent0pt}
\def\Note@HeadFont{\bf}
\def\Note@HeadToNumSep{}
\def\Note@HeadToSectionSep{\ }
\def\Note@SecOpenStyle{\bgroup(}
\def\Note@SecCloseStyle{)\egroup}
\def\Note@SecPunctuationStyle{.}
\def\Note@HeadParaSep{\hskip2pt}
\def\Note@AboveSpace{12pt plus2pt minus1pt}
\def\Note@BelowSpace{12pt plus2pt minus1pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Remark %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{remark}
\def\remark@ParaIndent{0pt}
\def\remark@StartIndent{\leftskip0pt}
\def\remark@Name{}
\def\remark@Font{\itshape}
\def\remark@HeadFont{\bf}
\def\remark@HeadToNumSep{}
\def\remark@HeadToSectionSep{\ }
\def\remark@SecOpenStyle{\bgroup}
\def\remark@SecCloseStyle{\egroup}
\def\remark@SecPunctuationStyle{}%{ \thesection .\global\advance\remarkcount by 1
%\the\remarkcount.}
\def\remark@HeadParaSep{\quad}
\def\remark@AboveSpace{8pt plus2pt minus2pt}
\def\remark@BelowSpace{8pt plus2pt minus2pt}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{quote}{\par\addvspace{6pt plus2pt}
\def\tabular{%\tablefont%
\def\@xarraycr{\@ifnextchar[\@argarraycr{\ifnum0=`{\fi}${}\cr\noalign{\vskip0pt}}}%
\let\@halignto\@empty\@tabular}%
\@namedef{tabular*}##1{%\tablefont%
\def\@xarraycr{\@ifnextchar[\@argarraycr{\ifnum0=`{\fi}${}\cr\noalign{\vskip0pt}}}%
\def\@halignto{to##1}\@tabular}%
      \quotefont\noindent\ignorespaces
}{%
      \par\addvspace{6pt plus2pt}}
%
\newenvironment{extract}{\par\addvspace{6pt plus2pt}
      \extractfont\noindent\ignorespaces
}{%
      \par\addvspace{6pt plus2pt}}
%
\def\vtfont{\fontsize{9}{10}\selectfont}
%
\usepackage{picinpar}
%
\newif\ifwindow
\global\windowfalse

\def\vt{\@ifnextchar[{\@vt}{\@@vt}}

\def\@vt[#1]{\global\windowtrue%
\vtfont\begin{window}[0,l,{#1},{}]
\noindent\ignorespaces}

\def\@@vt{\par\addvspace{12pt plus2pt minus2pt}%
\vtfont\noindent}

\def\endvt{\ifwindow\end{window}\fi\global\windowfalse
\par\addvspace{12pt plus2pt minus2pt}}

\def\glossname#1{\gdef\@glossnam{#1}}
\glossname{Glossary}
\def\glossheadfont{\fontsize{8}{10}\selectfont\itshape}
\def\glossdeffont{\fontsize{8}{10}\selectfont}
%
\NEWCmd\glosshead[1]{{\glossheadfont\noindent#1:}\ }
\NEWCmd\glossdefintion[1]{{\HANGin12pt\glossdeffont#1\par}}
%
\GlossaryENV{glossary}[1]{\par
\def\haone{#1}
\setlength{\parindent}{0pt}
\RENEWCmd\section{\@startsection {section}{1}{\z@}{-24pt \@plus -2pt \@minus -.2pt}{8pt}{\sectionfont}}%
\def\subsectionfont{\fontsize{8}{10}\selectfont\itshape}
\RENEWCmd\subsection{\@startsection{subsection}{2}{\z@}{-10pt \@plus -2pt \@minus -.2pt}{10pt}{\subsectionfont}}%
\ifx\haone\@empty
\else
\section*{#1}\fi
\reset@font}
{\par}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Special Definitions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newbox\tempboxa
\newdimen\calcwd
\newdimen\calcht
%
\def\<ss="#1" Place="#2" Matter="#3">{
\setbox\tempbox\vbox{#3}
\setbox\tempboxa\hbox{#3}
\calcwd\wd\tempboxa
\calcht\ht\tempbox
\advance\calcht2pt
\def\PlaceTOP{#3\llap{\raise\calcht\hbox to \calcwd{\hfill\vbox to 0pt{\hbox{$#1$}}\hfill}}}
\def\PlaceBOTTOM{#3\llap{\lower1pt\vbox to 0pt{\hbox to \calcwd{\hfill$#1$\hfill}}}}
\def\PlaceMIDDLE{#3\llap{\raise.5\calcht\vbox to 0pt{\hbox to \calcwd{\hfill$#1$\hfill}}}}
\def\Place{\@latex@error{^^J
Place must contain a Placement Pacifier
^^J}{}
}
\MakeCommand{Place#2}
}
%
\let\@@overline\overline
\def\overline#1{%
  \relax
  \ifmmode\@@overline{#1}%
  \else $\@@overline{\hbox{#1}}\m@th$\relax\fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USE PACKAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\usepackage{Bookmarks}
\usepackage{amssymb}
\usepackage{amsmath}

\def\eurl#1{{\insert\footins{\vskip1.5pt\noindent\footnotesize{{\itshape URL:}}\ #1.\vskip2.8pt}}}
%
\def\eurls#1{{\insert\footins{\vskip2pt\noindent\footnotesize{{\itshape URLs:}}\ #1.\vskip2.8pt}}}
%
\def\enddocument{%
\PDFTitle{\@forPDFTitle}{}
\immediate\write\bookmarks{\string\gdef\string\@lastpage{\the\c@page}}
   \@enddocumenthook
   \@checkend{document}
   \clearpage
\immediate\closeout\bookmarks
   \begingroup
     \if@filesw
       \immediate\write\@mainaux{\string\questionmark{}}
       \immediate\closeout\@mainaux
       \let\@setckpt\@gobbletwo
       \let\@newl@bel\@testdef
       \@tempswafalse
       \makeatletter \input\jobname.aux
     \fi
     \@dofilelist
     \ifdim \font@submax >\fontsubfuzz\relax
       \@font@warning{Size substitutions with differences\MessageBreak
                  up to \font@submax\space have occured.\@gobbletwo}%
     \fi
     \@defaultsubs
     \@refundefined
     \if@filesw
       \ifx \@multiplelabels \relax
         \if@tempswa
           \@latex@warning@no@line{Label(s) may have changed.
               Rerun to get cross-references right}%
         \fi
       \else
         \@multiplelabels
       \fi
     \fi
   \endgroup
   \deadcycles\z@\@@end}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\ifshowstats
\def\TRIMHeight{240mm}
\def\TRIMWidth{165mm}
\def\PAGEHeight{208mm}
\def\PAGEWidth{161mm or 36pc}
\def\COLUMNWidth{17pc}
\def\COLUMNSep{2pc}
\def\REFERENCEStyle{Numbered}
%
\def\ThreeSpaces{\space\space\space}
\def\FourSpaces{\space\space\space\space}
\def\FiveSpaces{\space\space\space\space\space}
\def\SixSpaces{\space\space\space\space\space\space}
\def\SevenSpaces{\space\space\space\space\space\space\space}
\def\EightSpaces{\space\space\space\space\space\space\space\space}
\def\TwentySpaces{\EightSpaces\EightSpaces\FourSpaces}

\read-1 to \tempskip
\fi
%
\vbadness=10000
\hbadness=10000
\vfuzz=10001pt
\hfuzz=10001pt
%
%%%%%%%%%%%%% Remove Package Warning for array.sty %%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand\Nnormalsize{%
   \@setfontsize\normalsize{11}{13}%\@xipt{13}
   \abovedisplayskip 7\p@ \@plus2\p@ \@minus.5\p@
   \abovedisplayshortskip \abovedisplayskip
   \belowdisplayshortskip \abovedisplayskip
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 29-01-05

\def\etal{{\it{ET~AL.}}}
%
\def\Tblequation{$\displaystyle}
\def\endTblequation{$\hfill \@eqnnum}
\def\unnumeq{$\displaystyle}
\def\endunnumeq{$}

\def\tabular{\tbzdbondtrue\tablefont%
\let\equation\Tblequation
\let\endequation\endTblequation
\let\[\unnumeq
\let\]\endunnumeq
\def\@xarraycr{\@ifnextchar[\@argarraycr{\ifnum0=`{\fi}${}\cr\noalign{\vskip2pt}}}%
\let\@halignto\@empty\@tabular}
\@namedef{tabular*}#1{\tbzdbondtrue\tablefont%
\let\equation\Tblequation
\let\endequation\endTblequation
\let\[\unnumeq
\let\]\endunnumeq
\def\@xarraycr{\@ifnextchar[\@argarraycr{\ifnum0=`{\fi}${}\cr\noalign{\vskip2pt}}}%
\def\@halignto{to#1}\@tabular}

\def\vt{\@ifnextchar[{\@vt}{\@@vt}}%

\def\@vt[#1]{\addvspace{12pt}\setbox\tempbox\hbox{#1}%
\vtfont
\vbox to 10pt{\box\tempbox}%
\vspace*{-12pt}%
\hangindent72.2pt\hangafter-9%
\vtfont\noindent\ignorespaces}%

\def\@@vt{\addvspace{12pt}%\setbox\tempbox\hbox{#1}%
\vtfont\noindent\ignorespaces}%

\def\endvt{\par}%
\def\query#1{\hyperlink{#1}{$^{\underline{\rm{{\blue{#1}}}}}$}}
\def\answer#1#2{\hypertarget{#1}{}{\noindent\blue #1: Author:~{#2\par}}\black}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START S2T CONVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\abscopyright#1{%
\gdef\@abscopyright{#1}}%
\abscopyright{}

\def\abscopyright{\abscoptrue{~\@jname}~{\ifearlyview\else\@jvol:~{\@firstpage--\@lastpage;~\@jyear}\fi\quad}\ifx\@checkcopyrightline\@empty\else$\protect\copyright$~\ifx\@cjyear\@empty\@jyear\else\@cjyear\fi~\@copyrightline\fi}%

\def\Frenchabscopyright{\abscoptrue{~\@Fjname}~{\ifearlyview\else\@jvol:~{\@firstpage--\@lastpage;~\@jyear}\fi\quad}\ifx\@checkFrenchcopyrightline\@empty\else\copyright~\ifx\@cjyear\@empty\@jyear\else\@cjyear\fi~\@Frenchcopyrightline\fi}%

\def\queryhead#1{\gdef\@queryhead{#1}}
\queryhead{Author Query Form (CJS/\@aid)}

\newenvironment{querypage}{\par%
                           \newpage\onecolumn\thispagestyle{empty}
                           {{\fontsize{15}{18}\selectfont\centerline{\bfseries \@queryhead}\par}\addvspace{24pt}}
                           {\bfseries\noindent Special Instruction: Author please include responses to queries with your other corrections and return by e-mail.\par\vspace*{12pt}}
                           }
               {\par}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hypothesis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NewEnvironment{hypothesis}
\def\hypothesis@ParaIndent{0pt}
\def\hypothesis@StartIndent{\leftskip0pt}
\def\hypothesis@Name{}
\def\hypothesis@Font{\itshape}
\def\hypothesis@HeadFont{\bf}
\def\hypothesis@HeadToNumSep{}
\def\hypothesis@HeadToSectionSep{\ }
\def\hypothesis@SecOpenStyle{\bgroup}
\def\hypothesis@SecCloseStyle{\egroup}
\def\hypothesis@SecPunctuationStyle{}%{ \global\advance\definitioncount by 1 \the\definitioncount.}
\def\hypothesis@HeadParaSep{\quad}
\def\hypothesis@AboveSpace{8pt plus2pt minus2pt}
\def\hypothesis@BelowSpace{8pt plus2pt minus2pt}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\lrborder#1#2#3{\smash{\raise#3pt\hbox{{\hskip#2pt}%
{\noindent\vrule height#1pt width.5pt depth0pt}%
{\hskip-#2pt}}}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@AFFIL{%
{\vtop{\hbox{{AUTHOR}\par}
\hbox{AUTHOR DEPT}\par
\hbox{ADD}}\par}}

\def\@HIST{{\fontsize{9}{12}\selectfont\itshape\noindent\ifx\@received\@empty\else Received~\@received\par\fi%
\noindent\ifx\@accepted\@empty\else Accepted~\@accepted\par\fi%
\noindent\ifx\@revised\@empty\else Revised~\@revised\par\fi}}




\def\CJShistory{%
\par\addvspace{-4.7pt}%
\setlength\parskip{0pt}
\noindent\vrule width\columnwidth height.5pt\par\addvspace{5pt}%
\def\rec##1##2##3{##1 ##2 ##3}
\def\acc##1##2##3{##1 ##2 ##3}
{\leftskip0pt\rightskip0pt plus1fill%
\smash{\hbox{\vtop{\@HIST}}}}\par\vspace*{40pt}}

\newenvironment{miscright}{\par\addvspace{-4.7pt}%
\setlength\parskip{0pt}
\noindent\vrule width\columnwidth height.5pt\par\addvspace{5pt}%
\leftskip0pt plus1fill\rightskip0pt%
\CJShistory
\def\dept##1{\fontsize{9}{12}\selectfont\itshape\par##1}%
\def\address##1{\fontsize{9}{12}\selectfont\itshape\par##1\par\addvspace{6pt}}%
\def\author##1{\fontsize{9}{12}\selectfont\par{\rm ##1}}%
\def\email##1{{{\selectfont$\!$:\,\,\,##1}}}%
}{\par\setlength\parskip{1pt}}



\newenvironment{pitmiscright}{\par\addvspace{11pt}
\leftskip0pt plus1fill\rightskip0pt%
\def\authorgroup##1{##1\par\addvspace{3pt}
\@DOInumber}
}{\par}








%%%%%%%%%%%%%%%%%%%%%%% QtoA for query link %%%%%%%%%%%%%%%%%%%%%%%%
\setlength\marginparwidth {11pt}
\setlength\marginparsep{5\p@}
\setlength\marginparpush{5\p@}
\def\red{\special{color cmyk 0 1 1 0}}

\newcommand{\leftmarginpar}[1]{%
\vadjust{\smash{\raise.2\baselineskip\vbox{\llap{\parbox[t]{\marginparwidth}{#1}\kern\marginparsep}}}}}
%
\newcommand{\rightmarginpar}[1]{%
\vadjust{\smash{\raise.2\baselineskip\vbox{\rlap{\parbox[t]{\marginparwidth}{\hbox to \columnwidth{\hfill\rlap{\hspace*{5pt}#1}}}\kern\marginparsep}}}}}
%
\newcommand{\Rightmarginpar}[1]{%
\vadjust{\smash{\raise.2\baselineskip\vbox{\rlap{\parbox[t]{\marginparwidth}{\hbox to \textwidth{\hfill\rlap{\hspace*{5pt}#1}}}\kern\marginparsep}}}}}
%
\newcommand{\mymarginpar}[1]{%
\vadjust{\smash{\raise.2\baselineskip\vbox{\llap{\parbox[t]{\marginparwidth}{#1}\kern\marginparsep}}}}}

\def\qtoa#1{\leftmarginpar{{\raggedright\fontsize{7}{9}\selectfont\ifdraft\smash{\special{color push}\black{ #1}\special{color pop}}\else\black\smash{ #1}\fi\endgraf}}}
\def\Qtoa#1{%
\rightmarginpar{{{\fontsize{7}{9}\selectfont\ifdraft\smash{\special{color push}\black{#1}\special{color pop}}\else\black\smash{#1}\fi}\endgraf}}%
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\hyperlink#1#2{\leavevmode\nobreak#2}%
\def\hypertarget#1{}
\def\mychar{}%
\def\url#1#2#3{#2}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usepackage{times}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END S2T CONVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\endinput
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                                               END OF CLASS FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
