using System;
using RDotNet.Internals;
using RDotNet;
using System.IO;
using Bio;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using RDotNet.NativeLibrary;
using System.Runtime.InteropServices;
namespace MitoDataAssembler
{
	/// <summary>
	/// R interface setup, creates an engine after setting the path approprioately, NOT thread safe.
    /// 
    /// This will be done by the PATH variable now, specifically R_HOME and R_LIB must be set
	/// </summary>
	public class RInterface
	{
		//TODO: This is dangerous, should be cleaned up once the configuration settings are better.
        internal const string R_LIB_ENV_DIR="R_LIB";
        internal const string R_HOME_ENV_DIR="R_HOME";
		private static RDotNet.REngine pEngine;
        private static object lockObject=new Object();
		public REngine CurrentEngine {
            get
            {
                if (pEngine == null) {
                    throw new Exception("Somehow tried to access engine without using constructor first");
                }
                return pEngine;
            }
		}
        //Hard coded dependency on R.NET 2.9
		public const string BROAD_RLIB_PATH="/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_2.9.2/lib64/R/lib";
		public const string UNIX_PATH = "PATH";
		public RInterface ()
		{
            //check for environmental variable to DLL file, note I custom recompiled this
            if (RInterface.pEngine == null)
            {
                lock (lockObject)
                {
                    //check for value once lock obtained, as may have been set already
                    if (RInterface.pEngine == null)
                    {
                        string dll = System.Environment.GetEnvironmentVariable(R_LIB_ENV_DIR);
                        if (String.IsNullOrEmpty(dll))
                        {
                            throw new Exception("R - library file was not set by environmental variable: " + R_LIB_ENV_DIR + ".\n Please set this variable to point to the directory with the library (libR.so, R.dylib or R.dll as needed).");
                        }
                        string r_home = System.Environment.GetEnvironmentVariable(R_HOME_ENV_DIR);
                        if (String.IsNullOrEmpty(r_home))
                        {
                            throw new Exception(R_HOME_ENV_DIR + " environmental variable is not set.  Please point this to your R Directory");
                        }
                        //change path 
                        System.Diagnostics.Debug.WriteLine(R_HOME_ENV_DIR + ": " + r_home);
                        System.Diagnostics.Debug.WriteLine(R_LIB_ENV_DIR + ":" + dll);

                     
                        if (System.Environment.OSVersion.Platform != PlatformID.Unix)
                        {
                            var envPath = Environment.GetEnvironmentVariable("PATH");
                            Environment.SetEnvironmentVariable("PATH", envPath + Path.PathSeparator + dll);
                            RInterface.pEngine = REngine.CreateInstance("RDotNet", "R");
                        }
                        else
                        {

                            RInterface.pEngine = REngine.CreateInstance("RDotNet", dll);

                        }

                        StartupParameter sp = new StartupParameter();
                        sp.Interactive = false;
                        sp.Slave = true;
                        sp.Verbose = false;
						sp.Quiet = true;
                        sp.SaveAction = StartupSaveAction.NoSave;
                        //THIS IS CRITICAL: See https://rdotnet.codeplex.com/workitem/70
                        var platform = Environment.OSVersion.Platform;
                        if (platform == PlatformID.Unix || platform == PlatformID.MacOSX)
                        {
                            System.Diagnostics.Debug.WriteLine("Removing R signal handlers");
                            IntPtr callBackPointer = RInterface.pEngine.DangerousGetHandle("R_SignalHandlers");
                            Marshal.WriteInt32(callBackPointer, 0);
                        }
                        RInterface.pEngine.Initialize(sp);
                    }
                }
            }
		}
        /// <summary>
        /// Make a PDF from the given x,y data.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="pdfName"></param>
        /// <param name="title"></param>
        /// <param name="xlab"></param>
        /// <param name="ylab"></param>
        /// <param name="AdditionalCommands"></param>
        /// <param name="?"></param>
        /// <param name="optArgs"></param>
        public void PlotPDF(IEnumerable<double> x, IEnumerable<double> y, string pdfName,  string title = null, string xlab = null, string ylab = null, List<string> AdditionalCommands = null,List<OptionalArgument> optArgs=null)
        {         
            int xl = x.Count();
            int yl = y.Count();
            if (xl != yl)
            {
                throw new Exception("Tried to plot vectors of unequal length!");
            }
            else if (xl == 0)
            {
                throw new Exception("Tried to plot zero length vector!");
            }
            List<OptionalArgument> args = new List<OptionalArgument>() { new OptionalArgument("main", title), new OptionalArgument("xlab", xlab), new OptionalArgument("ylab", ylab) };
                if (optArgs!=null)
            {
                args.AddRange(optArgs);
            }
            var engine = pEngine;
            lock (engine)
            {
                var yvals = engine.CreateNumericVector(y);
                var xvals = engine.CreateNumericVector(x);
                engine.SetSymbol("xvals", xvals);
                engine.SetSymbol("yvals", yvals);
                engine.Evaluate("library(grDevices)");
                string pdfCMD = "pdf(\"" + pdfName + "\")";
                engine.Evaluate(pdfCMD);
                var locOptArgs = args.Where(z => z.HasValue()).Select(b => b.ValuePair()).ToArray();
                string command = "plot(xvals,yvals";
                if (locOptArgs.Length > 0)
                {
                    command += "," + String.Join(",", locOptArgs);
                }
                command += ")";
                engine.Evaluate(command);
                if (AdditionalCommands != null)
                {
                    foreach (var cmd in AdditionalCommands)
                    {
                        engine.Evaluate(cmd);
                    }
                }
                engine.Evaluate("dev.off()");
            }            
        }

        public class OptionalArgument { 
            protected string name,value;
            public OptionalArgument(string name, string value) { this.name = name; this.value = value; }
            public bool HasValue() { return !String.IsNullOrEmpty(value); }
            public virtual string ValuePair() { return name + "=\"" + value + "\""; }        
        }

        public class NumericOptionalArgument : OptionalArgument
        {
            public NumericOptionalArgument(string s, string v) : base(s, v) { }
            public override string ValuePair()
            {
                return name + "=" + value + ""; 
            }
        }
        
        public void HistPDF(IEnumerable<double> x, string pdfName,int breaks=10, string title=null,string xlab=null,string ylab=null,List<string> AdditionalCommands=null,List<OptionalArgument> optArgs=null)
        {           
            List<OptionalArgument> args = new List<OptionalArgument>() { 
                                                new NumericOptionalArgument("breaks", breaks.ToString()), 
                                                new OptionalArgument("main", title), 
                                                new OptionalArgument("xlab", xlab), 
                                                new OptionalArgument("ylab", ylab) };
            if (optArgs!=null)
            {
                args.AddRange(optArgs);               
            }
            var engine = pEngine;
            lock (engine)
            {
                var xvals = engine.CreateNumericVector(x);
                engine.SetSymbol("xvals", xvals);
                engine.Evaluate("library(grDevices)");
                string pdfCMD = "pdf(\"" + pdfName + "\")";
                engine.Evaluate(pdfCMD);
                var locOptArgs = args.Where(z => z.HasValue()).Select(b => b.ValuePair()).ToList();
                string command = "hist(xvals,col=\"blue\"";
                if (locOptArgs.Count>0)
                {
                    command += "," + String.Join(",", locOptArgs);
                }    
                command += ")";
                engine.Evaluate(command);
                if (AdditionalCommands != null)
                {
                    foreach (var cmd in AdditionalCommands)
                    {
                        engine.Evaluate(cmd);
                    }
                }
                engine.Evaluate("dev.off()");
            }            
        }

		private static void setupMacOSX()
		{
            
			var rHome = "/Library/Frameworks/R.framework/Resources";
			Environment.SetEnvironmentVariable("R_HOME", rHome);
			string rLibDirec = "/Library/Frameworks/R.framework/Resources/lib";
			var oldPath = System.Environment.GetEnvironmentVariable(UNIX_PATH);
			if (Directory.Exists(rLibDirec) == false)
				throw new DirectoryNotFoundException(string.Format("Could not found the specified path to the directory containing R.dylib: {0}", BROAD_RLIB_PATH));
			var newPath = string.Format("{0}{1}{2}", rLibDirec, System.IO.Path.PathSeparator, oldPath);
			System.Environment.SetEnvironmentVariable("PATH", newPath);
			System.Environment.SetEnvironmentVariable ("R_HOME", rHome);
		}
		private static void SetupBroadR()
		{
			var oldPath = System.Environment.GetEnvironmentVariable(UNIX_PATH);
			if (Directory.Exists(BROAD_RLIB_PATH) == false)
				throw new DirectoryNotFoundException(string.Format("Could not found the specified path to the directory containing R.dll: {0}", BROAD_RLIB_PATH));
			var newPath = string.Format("{0}{1}{2}", BROAD_RLIB_PATH, System.IO.Path.PathSeparator, oldPath);
			System.Environment.SetEnvironmentVariable("PATH", newPath);
			System.Environment.SetEnvironmentVariable ("R_HOME", "/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_2.9.2/");
		}
		private static void setupWindowsR()
		{
			var oldPath = System.Environment.GetEnvironmentVariable(UNIX_PATH);
			var rPath = System.Environment.Is64BitProcess ? @"C:\Program Files\R\R-3.0.1\bin\x64" : @"C:\Program Files\R\R-3.0.1\bin\i386";
			if (Directory.Exists(rPath) == false)
				throw new DirectoryNotFoundException(string.Format("Could not found the specified path to the directory containing R.dll: {0}", rPath));
			var newPath = string.Format("{0}{1}{2}", rPath, System.IO.Path.PathSeparator, oldPath);
			System.Environment.SetEnvironmentVariable("PATH", newPath);
		}
        private static void setupPath()
        {
            if (RInterface.pEngine == null)
            {
                //var platform = Environment.OSVersion.Platform;
                var platform = Bio.CrossPlatform.Environment.GetRunningPlatform();
                switch (platform)
                {
                    case Bio.CrossPlatform.Environment.Platform.Windows:
                        setupWindowsR();
                        break; // R on Windows seems to have a way to deduce its R_HOME if its R.dll is in the PATH
                    case Bio.CrossPlatform.Environment.Platform.Mac:
                        setupMacOSX();
                        break;
                    case Bio.CrossPlatform.Environment.Platform.Linux:
                        SetupBroadR();
                        break;
                    default:
                        throw new NotSupportedException(platform.ToString());
                }
            }
        }

	}
}

