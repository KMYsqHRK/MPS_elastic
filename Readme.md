# MPS_elastic

This repository is being developed to simulate elastic bodies using the particle method. In the future, we aim to coupled with fluids.<br>
These codes were created using the following repositories as references.<br>
[MPS-basic](https://github.com/MPS-Basic/MPS-Basic.git)

## Requirement
### Execution
- Git
- cmake (newer than 3.9)
- C++ 17 compiler
- OpenMP 5.0 and above (optional)

### Dependencies
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

## Execution
### Build
1. Generate build system
	```bash
	cmake -S . -B build
	```
1. Execute build
	```bash
	cmake --build build
	```

### Execution
#### Windows
1. Create output directory if not exist
	```powershell
	New-Item -ItemType Directory -Path result/dambreak -Force
	```
1. Remove old output files if exist
	```powershell
	Remove-Item -Path $outputDir/* -Force -Recurse
	```
3. Run simulation
	```powershell
	./build/mps.exe --setting input/dambreak/settings.yml --output result/dambreak 2> result/dambreak/error.log | Tee-Object -FilePath "result/dambreak/console.log"
	```

#### Linux/Mac
1. Create output directory if not exist
	```bash
	mkdir -p result/dambreak/
	```
1. Remove old output files if exist
	```bash
	rm -rf result/dambreak/*
	```
1. Run simulation
	```bash
	./build/mps --setting input/dambreak/settings.yml --output result/dambreak 2> result/dambreak/error.log | tee result/dambreak/console.log
	```

### Checking Results
Result files will be written in `result/dambreak/vtu/***.vtu`.
Open these files in [ParaView](https://www.paraview.org/) to see the result.

