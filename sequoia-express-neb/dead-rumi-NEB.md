## Updating DEAD rumi version of SEQuoia Express Toolkit to accomodate NEB libraries 

### nextflow.config 
line 36: ADD neb = false

### src/neb.json
ADD file

### main.nf
1. line 73 ADD: if(!params.neb){ summary['Library Type'] = 'SEQuoia Express' } else { summary['Library Type'] = 'NEB' }
2. line 135 ADD: 
	if(!params.neb){			
		"""
		export RUST_LOG=info	
		dead -c /opt/biorad/src/2dcomplete.json -a DefaultParser -o ./ -i $sample_id $reads
		"""
	} else {
		"""
		export RUST_LOG=info	
		dead -c /opt/biorad/src/neb.json -a DefaultParser -o ./ -i $sample_id $reads
		"""
	}
 3. line 165 ADD: 	
 	if(!params.neb){
		cutter = "-u 1"
	} else {
		cutter = ""
	}
	if(params.noTrim){
		cutter = ""
	}
	
  4. line 188 ADD:
  	if (!params.neb && !params.noTrim){
			cutter = cutter+" -U 8"
	} else if(params.neb) {
			cutter = cutter+" -U 11"
	}


