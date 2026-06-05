# Python workflow

`tools/negaccel_workflow.py` owns case authoring and scan orchestration while the C++ simulator owns execution of one resolved case JSON.

## Preferred authoring flow

Edit the organized Python-side authoring file first, then materialize the canonical runtime JSON from it:

```bash
python3 tools/negaccel_workflow.py author-case negaccel-authoring.example.json
```

This writes `MTF_authoring_example/MTF_authoring_example.json` by default. The authoring model is grouped by geometry, boundary conditions, particle source, magnetic field, gas interactions, surface interactions, run controls, and outputs.

You can override authoring values directly from the CLI too:

```bash
python3 tools/negaccel_workflow.py author-case negaccel-authoring.example.json \
  --case-tag MTF_fast_test \
  --set run.particleCount=25000 \
  --set boundaryConditions.gridVoltagesVolts.extractionGrid=7000
```

## Edit and run one case

If you want to inspect or hand-edit the resolved runtime contract directly, start from the runtime example case, change the fields you care about, then run the simulator on that JSON:

```bash
cp negaccel-config.example.json my_case.json
```

Typical edits in `my_case.json` are:

- `simulation.particleCount` to change the Monte Carlo particle count
- `boundaryConditions.electrodes[*].voltageVolts` to change grid voltages
- `physics.stripping.minimumZMeters` to move the secondary-generation start plane
- `physics.surfaceCollisions.minimumImpactZMeters` to move the surface-collision generation cutoff
- `gasDensity.profiles[*].source.path` to point to a different tabulated density file

Then build and run:

```bash
source ./setup_environment.sh
make -j4 runtest_new_v2
./runtest_new_v2 my_case.json
```

`runtest_new_v2` accepts either a full JSON path or a case stem. It writes outputs next to the case JSON.

## Build one case

This command works on the canonical runtime JSON template directly:

```bash
python3 tools/negaccel_workflow.py build-case \
  negaccel-config.example.json \
  --case-tag MTF_custom_case \
  --set simulation.particleCount=25000 \
  --set boundaryConditions.electrodes[0].voltageVolts=7000
```

By default this writes `MTF_custom_case/MTF_custom_case.json`.

## Expand a scan

```bash
python3 tools/negaccel_workflow.py expand-scan negaccel-scan.example.json
```

This creates one folder per generated case and writes a `scan-manifest.json` under the configured scan output directory.

The scan template can point either to the canonical runtime JSON or to the higher-level authoring JSON. `negaccel-scan.example.json` now demonstrates the preferred authoring-first flow, so scan overrides and `scan.parameterPath` are applied to the organized authoring fields before each case is materialized into the runtime contract.

## Run a scan

```bash
python3 tools/negaccel_workflow.py run-scan \
  negaccel-scan.example.json \
  --simulator ./runtest_new_v2
```

The Python tool expands the scan first, then runs the simulator once per generated case JSON.
