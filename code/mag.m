(* Mathematica utility file for masters thesis *)
pltopts = {Frame -> True, FrameStyle -> Black, GridLines -> {Automatic, Automatic}, Axes -> False};
FPlot[args__] := Plot[args, Evaluate[pltopts]];
