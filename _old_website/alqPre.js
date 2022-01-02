function getParents(n){
	let set = genPre(n);
	set = sieveOversized(set, n)
	set = sieveNonFactors(set)
	set = sieveNonParent(set)
	return set.map(ele => ele.reduce((e,acc) => e*acc))
}
console.log(getParents(396));

//Generates 
function genPre(n){
	let lNums = []
	for(let i = 2; i <= n/2; i++) lNums.push(i);
	
	let ps = bwPowerSet(lNums);
	ps.forEach(ele => ele.unshift(1))
	return ps
	
}

function sieveOversized(set, n){
	let precs = []

	for(let i = 0; i < set.length; i++){
		let sSum = set[i].reduce((a,b) => a+b, 0)
		let diff = n-sSum
		if(diff > n/2 && diff != n ){
			let temp = set[i].slice()
			temp.push(diff)
			precs.push(temp)
		}else if(diff==0){
			let temp = set[i].slice()
			precs.push(temp)
		}
	}
	return precs
}

function sieveNonFactors(set){
	let nSet =[]
	for(let i=0; i < set.length; i++){
		let e = set[i]
		let flag = true
		for(let j = e.length-1; j>=0; j--){
			
			let facs = factor(e[j])
			
			if(!facs.every(ele => e.includes(ele))){
				flag = false
			}
		}
		if(flag){
			nSet.push(e)
		}
	}
	return nSet
}

function sieveNonParent(set){
	let nSet =[]
	
	for(let i=0; i < set.length; i++){
		let e = set[i]
		
		let flag = true
		let orig = e[1] * e[e.length-1]
	
		for(let h = 1; h < e.length/2; h++){
			//console.log(e[1+h]* e[e.length-1-h])
			if(e[1+h]* e[e.length-h-1] != orig){
				
			}
		}
		

		if(flag && arraysEqual(factor(orig), e))nSet.push(e)
	}
	
	return nSet
}


https://github.com/trekhleb/javascript-algorithms/blob/master/src/algorithms/sets/power-set/bwPowerSet.js
function bwPowerSet(originalSet) {
  const subSets = [];

  // We will have 2^n possible combinations (where n is a length of original set).
  // It is because for every element of original set we will decide whether to include
  // it or not (2 options for each set element).
  const numberOfCombinations = 2 ** originalSet.length;

  // Each number in binary representation in a range from 0 to 2^n does exactly what we need:
  // it shows by its bits (0 or 1) whether to include related element from the set or not.
  // For example, for the set {1, 2, 3} the binary number of 0b010 would mean that we need to
  // include only "2" to the current set.
  for (let combinationIndex = 0; combinationIndex < numberOfCombinations; combinationIndex += 1) {
    const subSet = [];

    for (let setElementIndex = 0; setElementIndex < originalSet.length; setElementIndex += 1) {
      // Decide whether we need to include current element into the subset or not.
      if (combinationIndex & (1 << setElementIndex)) {
        subSet.push(originalSet[setElementIndex]);
      }
    }

    // Add current subset to the list of all subsets.
    subSets.push(subSet);
  }

  return subSets;
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