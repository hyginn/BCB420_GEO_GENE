# BCB420_GEO_GENE

Please use `./R/main.R`, and edit the three parameters at the top of the R script.

```
# json filename
inputFile <- "sample.json"

# output filename for RData
outputFile <- "sampletable.RData"

# log file
logFile <- "sample.log"
```

## JSON Format

Each experiment (comparison) is considered to be a JSON object, with the parameters "geoSeries", "description", "controlSamples", and "experimentalSamples". The values for each parameter should be a string. "geoSeries" describes the GEO Series identifier. "description" should be a unique description or descriptor for the experiment. "controlSamples" defines the normal sample(s), and "experimentalSamples" defines the experimental sample(s).

Sample JSON input:
```
[{
		"geoSeries": "GSE35330",
		"description": "Experiment 1",
		"controlSamples": ["GSM866121", "GSM866122", "GSM866124", "GSM866126"],
		"experimentalSamples": ["GSM866111", "GSM866116", "GSM866118", "GSM866119"]
	}
]
```
