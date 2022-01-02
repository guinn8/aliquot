//Yang & Pomerance algorithm for computing untouchable Alq. numbers

//Contains the sum of divisors for every odd number (1, x]
//sigmaM[i] = sigma((2i)+3)
//sigmaM[(i-3)/2] = sigma(i)
let sigmaM = [] 

//Represents all even numbers [1, x] 
//After processing characFunc[i] = 0 indicates that (i+1)*2 is touchable
//characFunc[i] = 1 indicates that (i+1)*2 is untouchable
//characFunc[(i/2)-1] = 0 indicates that i is touchable


let knowUntouch = [2,5,52,88,96,120,124,146,162,188,206,210,216,238,
 246,248,262,268,276,288,290,292,304,306,322,324,
 326,336,342,372,406,408,426,430,448,472,474,498,
 516,518,520,530,540,552,556,562,576,584,612,624,
 626,628,658]
 
 console.log(arr_diff(knowUntouch, driver(658)))



function driver(x){
	let characFunc = []
	//Calculates sigma(m) for all odd numbers (1, X]
	for(let m = 3; m <=x; m+=2){
		sigmaM[(m-3)/2] = sigma(m)
	}
	
	for(let i = 2; i<= x; i+=2){
		characFunc[(i/2)-1] = 1
	}
	//For each odd number m repersented in sigmaM this computes the recurrance
	//Base: t = 3 * sigma(m) -2m
	//Non-base: t = 2t + sigma(m) <<continue till t exceeds X>>
	//Each value of t is an even touchable number and is recorded in characFunc[t] = 0
	for(let i = 0; i < sigmaM.length; i++){
		let sigma_m = sigmaM[i] //Lookup of sigma(m)
		let m = (2*i)+3 //Recovers the odd number m such that sigma(m) is saved at index i 
		
		//If sigma_m is odd we can ignore that value as we are only ennumerating even values
		//Fun fact: sigma_m is odd for an even value of m iff m is a square or twice a square 
		if(sigma_m % 2 == 0){
			let t = 3*(sigma_m) - 2*m //init value of recurrance
			//Running the recurrance
			while(t <= x){
				
				characFunc[(t/2)-1] =0 
				t = 2 * t + sigma_m
			}
		
		}
		
		//This detects if m is a prime as the sigma(p) for any prime p is equal to p +1
		if(sigma_m == m+1){
			
			//The value of m+1 is recorded as a touched number because we know s(m^2) = m+1
			//When you square a prime p we know that numbers only proper factors
			characFunc[((m+1)/2)-1] = 0
		
		//Detect non-prime odd numbers m
		//s(m^2) is even
		}else if(m < Math.pow(x, 2/3)) {
			
			let mSq = Math.pow((m), 2)
			let s_mSq = sigma(mSq)- mSq
		
			if(s_mSq <= x){
				characFunc[(s_mSq/2)-1] = 0
			}
		}
	}
	
	return recUntouch(characFunc)
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

function recUntouch(characFunc){
	let untouched = [5]
	for(let i = 0; i < characFunc.length; i++){
		if(characFunc[i] == 1 ) {
			untouched.push((i+1)*2)
		}
	}
	return untouched
	
}

function arr_diff (a1, a2) {

    var a = [], diff = [];

    for (var i = 0; i < a1.length; i++) {
        a[a1[i]] = true;
    }

    for (var i = 0; i < a2.length; i++) {
        if (a[a2[i]]) {
            delete a[a2[i]];
        } else {
            a[a2[i]] = true;
        }
    }

    for (var k in a) {
        diff.push(k);
    }

    return diff;
}