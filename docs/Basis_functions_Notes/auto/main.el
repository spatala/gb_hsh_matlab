(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("wlscirep" "fleqn" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("ulem" "normalem") ("matlab-prettifier" "numbered" "framed")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "Introduction"
    "wlscirep"
    "wlscirep11"
    "lineno"
    "hyperref"
    "ulem"
    "xspace"
    "multicol"
    "setspace"
    "color"
    "mathtools"
    "bm"
    "tensor"
    "matlab-prettifier"
    "sidecap"
    "esvect")
   (TeX-add-symbols
    '("norm" 1)
    '("TTm" 1)
    '("CBred" 1)
    '("RBaseL" 1)
    '("Vect" 1)
    '("BaseL" 1)
    '("LTBasism" 2)
    '("LBasism" 2)
    '("LBasis" 2)
    '("VecSpaceT" 1)
    '("VecSpace" 1)
    '("RotMat" 3)
    '("Cred" 1)
    "Tref"
    "beginsupplement"
    "GrimmP"
    "BaseG"
    "colg"
    "colb"
    "colr"
    "ph")
   (LaTeX-add-bibliographies)
   (LaTeX-add-mathtools-DeclarePairedDelimiters
    '("ceil" "")
    '("floor" "")))
 :latex)

