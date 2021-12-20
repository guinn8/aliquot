function distinctPrimes(n){
    if(n % 2 != 0) n--

    for(p = n; p > 1; p--){
        if(prime(p)){
            for(q = 2; q < p; q++){
                if(prime(q)){
                    if(p+q==n) return [p, q]
                }
            }
        }
    }
}


function prime(n){
	let max = Math.floor(Math.sqrt(n))
	for(i = 2; i <= max; i++){
		if(n % i == 0){
			return false
		}
	}
	return true	
}