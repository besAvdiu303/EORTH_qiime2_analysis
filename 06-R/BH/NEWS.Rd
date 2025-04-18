\name{NEWS}
\title{News for Package \pkg{BH}}
\newcommand{\boost}{\href{https://www.boost.org}{\pkg{Boost}}}
\newcommand{\boostPkg}{\href{https://www.boost.org/doc/libs/release/libs/#1/}{\pkg{#1}}}
\newcommand{\boostPkgTwo}{\href{https://www.boost.org/doc/libs/release/libs/#1}{\pkg{#2}}}
\newcommand{\rfticket}{\href{https://r-forge.r-project.org/tracker/index.php?func=detail&aid=#1}{#1}}
\newcommand{\ghit}{\href{https://github.com/eddelbuettel/bh/issues/#1}{##1}}
\newcommand{\ghpr}{\href{https://github.com/eddelbuettel/bh/pull/#1}{##1}}

\section{Changes in version 1.87.0-1 (2024-12-17)}{
  \itemize{
    \item Upgrade to \boost 1.87.0, patched as usual to comment-out
    diagnostic suppression messages per the request of CRAN
    \item Switched to Authors@R
  }
}

\section{Changes in version 1.84.0-0 (2024-01-09)}{
  \itemize{
    \item Upgrade to \boost 1.84.0 
    \item Added \boostPkg{qvm} (requested in \ghit{97})
  }
}

\section{Changes in version 1.81.0-1 (2023-01-17)}{
  \itemize{
    \item Upgrade to \boost 1.81.0 (\ghpr{87})
    \item Added \boostPkg{url} (new in 1.81.0)
    \item Converted remaining \code{sprintf} to \code{snprintf}
    (\ghpr{90} fixing \ghit{89})
    \item Comment-out gcc warning messages in three files
  }
}
    
\section{Changes in version 1.78.0-0 (2021-12-14)}{
  \itemize{
    \item Incremental updates to README.md and CI setup
    \item Upgraded to \boost 1.78.0 (\ghpr{82} closing \ghit{80})
    \item Added \boostPkg{process} (closing \ghit{78})
    \item Added \boostPkg{lambda2} (closing \ghit{81})
  }
}

\section{Changes in version 1.75.0-0 (2020-12-12)}{
  \itemize{
    \item Removed file \code{NAMESPACE} as the package has neither R code,
    nor a shared library to load
    \item The file LICENSE_1_0.txt is now included (as requested in \ghit{73})
    \item Added new \boostPkg{beast} library (as requested in \ghit{74})
    \item Upgraded to \boost 1.75.0 (\ghpr{75})
  }
}

\section{Changes in version 1.72.0-3 (2019-12-28)}{
  \itemize{
    \item Remove old Boost installation when updating via script
  }
}

\section{Changes in version 1.72.0-2 (2019-12-16)}{
  \itemize{
    \item Bugfix release correction numeric/odeint/stepper/... renaming
  }
}

\section{Changes in version 1.72.0-1 (2019-12-15)}{
  \itemize{
    \item Upgraded to \boost 1.72.0 (plus the few local tweaks) (\ghpr{65})
    \item Applied the standard minimal patch with required changes, as
    well as the newer changeset for diagnostics pragma suppression.
    \item No longer install \boostPkg{filesystem} _explicitly_ though
    some files are carried in (\ghit{55})
    \item Added \boostPkg{mp11} (as requested in \ghit{62})
    \item Added \boostPkg{polygon} (as requested in \ghit{63})
  }
}

\section{Changes in version 1.69.0-1 (2019-01-07)}{
  \itemize{
    \item Upgraded to \boost 1.69.0 (plus the few local tweaks)
    \item Applied the standard minimal patch with required changes, as
    well as the newer changeset for diagnostics pragma suppression.
    \item Following a pre-release in December, maintainers of three
    packages affected by the 1.66 to 1.69 were contacted, and changes
    were prepared.
  }
}

\section{Changes in version 1.66.0-1 (2018-02-12)}{
  \itemize{
    \item Upgraded to \boost 1.66.0 (plus the few local tweaks)
    \item Added Boost \boostPkg{compute} (as requested in \ghit{16})
  }
}
    
\section{Changes in version 1.65.0-1 (2017-08-24)}{
  \itemize{
    \item Upgraded to \boost 1.64 and then 1.65 installed directly from
    upstream source with several minor tweaks (as before)
    \item Fourth tweak corrects a misplaced curly brace (see the
    Boost ublas GitHub repo and its issue #40)
    \item Added Boost \boostPkg{align} (as requested in \ghit{32})
    \item Added Boost \boostPkg{sort} (as requested in \ghit{35})
    \item Added Boost \boostPkg{multiprecision} by fixing a script typo
    (as requested in \ghit{42})
    \item Updated Travis CI support via newer \code{run.sh}
  }
}

\section{Changes in version 1.62.0-1 (2016-11-15)}{
  \itemize{
    \item Upgraded to \boost 1.62 installed directly from upstream source
    \item Added Boost \boostPkg{property_tree} as requested in \ghit{29}
    by Aydin Demircioglu
    \item Added Boost \boostPkg{scope_exit} as requested in \ghit{30} by
    Kirill Mueller
    \item Added Boost \boostPkg{atomic} which we had informally added
    since 1.58.0
  }
}

\section{Changes in version 1.60.0-2 (2016-05-06)}{
  \itemize{
    \item Added Boost \boostPkg{bimap} via GH pull request \ghpr{24} by
    Jim Hester
    \item Added Boost \boostPkg{icl} via GH pull request \ghpr{27} by
    Jay Hesselbert
    \item Added Boost \boostPkg{flyweight} as requested in GH ticket \ghit{26}
  }
}

\section{Changes in version 1.60.0-1 (2015-12-24)}{
  \itemize{
    \item Upgraded to \boost 1.60 installed directly from upstream source
    \item Added Boost \boostPkg{phoenix} as discussed in GH ticket \ghit{19}
  }
}

\section{Changes in version 1.58.0-1 (2015-05-21)}{
  \itemize{
    \item Upgraded to \boost 1.58 installed directly from upstream source
    \item Added Boost \boostPkg{multiprecision} as requested in GH ticket \ghit{12} based
      on rcpp-devel request by Jordi Molins Coronado
  }
}

\section{Changes in version 1.55.0-3 (2015-01-04)}{
  \itemize{
    \item Added Boost \boostPkg{graph} requested in GH ticket
    \ghit{9} by Dirk for \CRANpkg{RcppStreams}
  }
}

\section{Changes in version 1.55.0-2 (2015-01-03)}{
  \itemize{
    \item Added Boost \boostPkg{fusion} requested in GH ticket
    \ghit{7} by Dirk for \CRANpkg{RcppStreams}
  }
}

\section{Changes in version 1.55.0-1 (2014-12-21)}{
  \itemize{
    \item Upgraded to \boost 1.55 once again from the  
    \href{https://packages.debian.org/source/sid/boost1.55}{Debian tarball}
    \item Added Boost \boostPkg{geometry} requested in GH ticket \ghit{5} by 'teramonagi'
  }
}

\section{Changes in version 1.54.0-5 (2014-11-09)}{
  \itemize{
    \item Added Boost \boostPkg{circular_buffer} requested by Ben Goodrich for \CRANpkg{rstan}
  }
}

\section{Changes in version 1.54.0-4 (2014-08-29)}{
  \itemize{
    \item Added Boost \boostPkg{any} requested by Greg Jeffries for his \CRANpkg{nabor} package 
  }
}

\section{Changes in version 1.54.0-3 (2014-08-03)}{
  \itemize{
    \item Added Boost \boostPkg{heap} library which will be needed by the next
    version of \CRANpkg{RcppMLPACK}
  }
}

\section{Changes in version 1.54.0-2 (2014-04-14)}{
  \itemize{
    \item Bug fix to \code{lexical_cast.hpp} which now uses the test for
    INT128 which the rest of Boost uses, consistent with Boost 1.55 too.
  }
}

\section{Changes in version 1.54.0-1 (2014-04-07)}{
  \itemize{
    \item Upgraded to Boost 1.54.0
    \item Adjust build script \code{local/script/CreateBoost.sh} accordingly
    \item Renamed \code{generation_runge_kutta_cash_karp54_classic.hpp}
    to \code{generation_runge_kutta_cash_karp54_cl.hpp} to remain within
    100-character limit for \code{tar}; also adjust corresponding
    \code{#include} in \code{numeric/odeint/stepper/generation.hpp}
  }
}

\section{Changes in version 1.51.0-4 (2014-01-01)}{
  \itemize{
    \item Rewritten main package creation script to no longer scan for
    what \CRANpkg{bigmemory} and \CRANpkg{synchronicity} use, but rather
    explicitly copy over an (equivalent) list of explicitly-enumerated
    \boost libraries
    \item Repository moved from R-Forge to GitHub, scripts and layout
    adjusted accordingly
    \item Besides the implicitly expanded coverage by including the
    complete libraries, we also expanded from
    \code{math/doc/html/dist.html} to all of \boostPkg{math}.
  }
}

\section{Changes in version 1.51.0-3 (2013-10-19)}{
  \itemize{
    \item Expanded \boost coverage by adding \boost packages
    \boostPkg{spirit}, \boostPkg{foreach}, \boostPkg{algorithm},
    \boostPkg{iostreams}. 
  }
}

\section{Changes in version 1.51.0-2 (2013-07-29)}{
  \itemize{
    \item Expanded \boost coverage by adding \boost packages
    \boostPkg{random} and \boostPkg{unordered}. 
  }
}

\section{Changes in [unreleased] version 1.51.0-1 (2013-05-04)}{
  \itemize{
    \item Expanded \boost coverage by adding \boost packages \boostPkg{foreach},
    \boostPkgTwo{math/doc/html/dist.html}{statistical
      distributions}, and \boostPkg{iostreams}. This addresses,
    respectively, R-Forge issue tickets \rfticket{2527}, \rfticket{2533}
    and \rfticket{2768}. 
    \item The package can now also be used to compile CRAN package
    \CRANpkg{ExactNumCI} and \CRANpkg{rcppbugs}.
  }
}

\section{Changes in version 1.51.0-0 (2013-01-31)}{
  \itemize{
    \item First CRAN upload
    \item The package can be used to build \CRANpkg{bigmemory} and
    \CRANpkg{RcppBDT} as it already contains significant portions of Boost libraries
    algorithm, bind, concept, config, container, date_time, detail,
    exceptions, filesystem, functional, interprocess, intrusive, iterator,
    math, mpl, numeric, preprocessor, random, range, smart_ptr,
    type_traits, unordered.
  }
}

