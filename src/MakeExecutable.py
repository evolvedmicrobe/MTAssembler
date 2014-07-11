import os
# This command may need to run in the terminal ahead of time
# You may also need to add -liconv to the CC command and run it in the terminal if you 
# see something along the lines of:
# undefined reference to `libiconv_open'
os.system("setenv PKG_CONFIG_PATH $HOME/monoInstall/lib/pkgconfig/")
t=[x for x in os.listdir(os.getcwd()) if x.endswith((".dll",".exe"))]
cmd="mkbundle --static  --deps -o mtanalysis "+" ".join(t)
#above didn't work, order seems to matter
#cmd="mkbundle --deps -o curvefitter --static CurveFitterMonoGUI.exe  ShoOptimizer.dll Microsoft.Solver.Foundation.dll MatrixInterf.dll ShoArray.dll alglibnet2.dll MatrixArrayPlot.dll ZedGraph.dll GrowthCurveLibrary.dll "
print cmd
os.system(cmd)

