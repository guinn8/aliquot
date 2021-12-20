function factorPair(n){
	 
	let max = Math.floor(Math.sqrt(n))
	let factors = [] 
	for(i = 2; i <= max; i++){
		if(n % i == 0){
			factors.push([i, n/i])
		}
	}
	return factors
}

function factor(n){
	
	let max = Math.floor(Math.sqrt(n))
	let factors = [1] 
	for(i = 2; i <= max; i++){
		if(n % i == 0){
			factors.push(i)
			factors.push(n/i)
		}
	}
	return factors
	
}

function arraysEqual(a, b) {
	a = a.sort()
	b = b.sort()
  if (a === b) return true;
  if (a == null || b == null) return false;
  if (a.length != b.length) return false;

  // If you don't care about the order of the elements inside
  // the array, you should sort both arrays here.
  // Please note that calling sort on an array will modify that array.
  // you might want to clone your array first.

  for (var i = 0; i < a.length; ++i) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

function getOneParent(n, lim){
	let pathSet = [[]]
	let singletonPaths = []
	let hits = {}
	for(let i =1;i<=n ; i++){
	
		let val = i
		hits[val] = 1
		let path = [val]
		let trigger = 0
		while(val > 1 && !path.slice(0, path.length-1).includes(val) && trigger < lim){
			val = sumDiv(path[path.length-1])
			path.push(val)
			
			if(typeof hits[val] != "number" ) hits[val] = 1
			else{hits[val] += 1}
			trigger++
		}
		if(trigger == lim) path.push("maxed")
		pathSet.push(path)
	}
	for(let i = 0; i<pathSet.length; i++){
		for(let j = 0; j < pathSet[i].length; j++){
			if(hits[pathSet[i][j]] == 1){
				singletonPaths.push([pathSet[i], pathSet[i][j]])
				break
			}
		}
	}
	
	return singletonPaths
}

function sumDiv(m){
	let max = Math.sqrt(m)
	let acc = 1 
	for(i = 2; i <= max; i++){
		if(m % i == 0){
			if(i != (m/i)) acc = acc + i + (m/i)
			else acc = acc + i
		}
	}
	return acc
}