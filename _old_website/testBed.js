let AllT_a = {}
let map = {}

//Pom/Pol set construction
function T_a(y, x){
	let A_y = math.lcm.apply(null, range(1,y+1))//

	let factors = factor(A_y)

	for(let i = 0; i < factors.length; i++){
		let a = factors[i]
		let n = 1
		let T_a = []
		let testN
		do{
			testN = math.gcd(n, A_y)
			if(testN == a) T_a.push(n)
			n++
		}while(n <= x)
		AllT_a[a] = T_a
	}
	return AllT_a
}

function checkMap(set){
	for(const ta in set){
		for(const val in set[ta]){
			if (typeof map[set[ta][val]] !== 'undefined') map[set[ta][val]].push(ta)
			else  map[set[ta][val]] = [ta]
			//console.log("ta: "+ ta + " val: " + set[ta][val])
		}

	}
	return map
}
//Test of claim that s maps Ta to Ta forall a that divides Ay. Asymp. true if n > e^e^^y as y -> infny
//The plan is iterate s(n) for everything in T_a and see what percentage of the time it lands back in T_a

//It seems to be true that a growing percentange of n in T_a is mapped back to T_a with greater values of x
//This effect does appear to reach some limiting value for the limited test cases
//The effect of changing the value of y appears to be much more dramatic in increasing the 
//average amount mapped back
function claim1(T_a){
	let acc = 0 
	for(let i = 0; i < T_a.length; i++){
		
		let iter = s(T_a[i])
		if(T_a.includes(iter)) acc++
	}
	return (acc/T_a.length)
}

//Test that sigma(n)/n is approx. equal to sigma(a)/a for n in T_a
//Asymp. true as y to infinity

//This also certainly seesm to be true. The values approach quite closly for even y = 20
function claim2(a, T_a){
	let aFrac = sigma(a)/a
	let fracVals = []
	for(let i = 0; i < T_a.length; i++){
		let n = T_a[i]
		let testFrac = sigma(n)/n
		fracVals.push(testFrac)
	}
	return [aFrac, fracVals]
}
      



//Riele lemma 9.2
//Im just doing this for sigma
function lemma9_2(d, x){
	let acc = 0
	for(let n = 1; n <= x; n++){
		let sig_n = sigma(n)
		console.log("n"+ n + " "+sig_n)
		if(sig_n % d != 0){
			acc++
		}
	}
	console.log("The number of positive integers n st d does not divide sigma(n) is " + acc)
	return acc
	
}

const range = (start, end) => {
    const length = end - start;
    return Array.from({ length }, (_, i) => start + i);
}

function sigma(m){
	let max = Math.sqrt(m)
	let acc = 1 
	for(i = 2; i <= max; i++){
		if(m % i == 0){
			if(i != (m/i)) acc = acc + i + (m/i)
			else acc = acc + i
		}
	}
	return acc + m
}

function factor(m){
	let max = Math.sqrt(m)
	let factors = [1, m] 
	for(i = 2; i <= max; i++){
		if(m % i == 0){
			if(i != (m/i)){
				factors.push(i)
				factors.push(m/i)
			}
			else factors.push(i)
		}
	}
	factors = factors.sort(function(a, b){return a-b})
	return factors
}

function s(m){
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