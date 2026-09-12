# Runtime requirements

JLinAlg requires JDK 17 or newer for source builds. The Gradle wrapper downloads
its pinned dependencies on the first build, including JDistlib 0.10.2, and
verifies the JDistlib JAR's SHA-256 digest. The executable JAR bundles its Java
dependencies; it does not install vendor GPU drivers or runtime libraries.

CPU-only use needs no GPU software. Optional accelerated backends require
separately installed, compatible native components:

- NVIDIA CUDA: [CUDA Toolkit 12.6](https://developer.nvidia.com/cuda-12-6-0-download-archive),
  a compatible NVIDIA driver, and the CUDA runtime, cuBLAS, and NVRTC libraries.
  On Windows, make the Toolkit's runtime DLL directory discoverable on `PATH`.
  CUDA 13 DLLs do not substitute for the required CUDA 12 libraries.
- OpenCL: an FP64-capable implementation supplied by the hardware vendor.
  See the [Khronos OpenCL page](https://www.khronos.org/opencl/).
- Vulkan: an FP64-capable Vulkan driver for the device.
  See the [LunarG Vulkan page](https://vulkan.lunarg.com/sdk/home).

For commands supporting backend selection, use `--backend cpu`, `--backend cuda`,
`--backend opencl`, or `--backend vulkan`. An explicitly requested GPU backend
fails if unavailable. Availability and benefit depend on the operation and
matrix size; see the [compute vignette](vignettes/formulas-and-backends.md).
