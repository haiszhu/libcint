class OpenblasSinglethread < Formula
  desc "Personal single-threaded OpenBLAS build"
  homepage "https://www.openblas.net/"
  url "https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v0.3.29.tar.gz"
  sha256 "38240eee1b29e2bde47ebb5d61160207dc68668a54cac62c076bb5032013b1eb"
  # to install: brew install --build-from-source --formula --verbose ./openblas-singlethread.rb
  # lsof /opt/homebrew/var/homebrew/locks/m4.formula.lock
  # kill -9 <PID>
  # export LDFLAGS="-L/opt/homebrew/opt/openblas-singlethread/lib"
  # export CPPFLAGS="-I/opt/homebrew/opt/openblas-singlethread/include"
  # gcc your_app.c -I/opt/homebrew/opt/openblas-singlethread/include -L/opt/homebrew/opt/openblas-singlethread/lib -lopenblas -o your_app
  # The main license is BSD-3-Clause. Additionally,
  # 1. OpenBLAS is based on GotoBLAS2 so some code is under original BSD-2-Clause-Views
  # 2. lapack-netlib/ is a bundled LAPACK so it is BSD-3-Clause-Open-MPI
  # 3. interface/{gemmt.c,sbgemmt.c} is BSD-2-Clause
  # 4. relapack/ is MIT but license is omitted as it is not enabled
  license all_of: ["BSD-3-Clause", "BSD-2-Clause-Views", "BSD-3-Clause-Open-MPI", "BSD-2-Clause"]
  head "https://github.com/OpenMathLib/OpenBLAS.git", branch: "develop"

  # Don't link into HOMEBREW_PREFIX; use manually via LDFLAGS/CPPFLAGS
  keg_only "personal single-threaded OpenBLAS"

  depends_on "gcc" # for gfortran
  fails_with :clang

  def install
    ENV.runtime_cpu_detection
    ENV.deparallelize # build is parallel by default, but setting -j confuses it
    # The build log has many warnings of macOS build version mismatches.
    ENV["MACOSX_DEPLOYMENT_TARGET"] = MacOS.version.to_s if OS.mac?
    ENV["DYNAMIC_ARCH"] = "1"
    ENV["USE_OPENMP"]   = "0"
    ENV["USE_THREAD"]   = "0"
    ENV["NUM_THREADS"]  = "1"
    # See available targets in TargetList.txt
    ENV["TARGET"] = case Hardware.oldest_cpu
                    when :arm_vortex_tempest then "VORTEX"
                    when :westmere          then "NEHALEM"
                    else Hardware.oldest_cpu.upcase.to_s
                    end
    # Apple Silicon does not support SVE
    # https://github.com/xianyi/OpenBLAS/issues/4212
    ENV["NO_SVE"] = "1" if Hardware::CPU.arm?

    opts = "-O3 -march=native -funroll-loops -fomit-frame-pointer"
    ENV["NO_SHARED"] = "1"
    
    # Version number issue on linux...
    ENV["VERSION"] = ""

    # Must call in two steps, no "shared" target
    system "make", "VERSION=", "CC=#{ENV.cc}", "FC=gfortran", "COMMON_OPT=#{opts}", "FCOMMON_OPT=#{opts}", "libs", "netlib"
    system "make", "PREFIX=#{prefix}", "install"

  end

  test do
    (testpath/"test.c").write <<~EOS
      #include <cblas.h>
      int main() {
        double x[1] = {1}, y[1] = {1};
        cblas_ddot(1, x, 1, y, 1);
        return 0;
      }
    EOS
    system ENV.cc, "test.c", "-I#{include}", "-L#{lib}", "-lopenblas", "-o", "test"
    system "./test"
  end

  def caveats
    <<~EOS
      To compile and link your own code against this single‑threaded OpenBLAS:

        gcc your_app.c \\
          -I#{opt_include} \\
          -L#{opt_lib} -lopenblas \\
          -o your_app

      (#{opt_include} expands to $(brew --prefix openblas-singlethread)/include)
    EOS
  end
  
end
