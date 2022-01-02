

SQ = 1280*2
SIZE = 80
OFF = 40
let slider;
let primes = []
let button;
let seq = [];

function setup() {
    createCanvas(SQ, SQ);
	slider = createSlider(0, 255, 10, 1);
	slider.position(SQ/2, SQ+10);
	slider.style('width', '80px');
	slider.hide()
	
	
	input = createInput();
	input.position(SQ/2, SQ+10);
	input.hide()
  
	button = createButton('submit');
	button.position(SQ/2, SQ+50);
	button.mousePressed(runSeq);
	button.hide()
	
	
	
	
	radio = createRadio();
	radio.option('Show proper factors',1);
	radio.option('other',2);
	radio.style('width', '150px');
	textAlign(CENTER);
	

	l = 2
	while(primes.length < SIZE){
		
		if(prime(l)){
			primes.push(l)
		} 
		l= l+1
		
	}
	
}



function draw() {
	clear();
	
	if(radio.value() == 1){
		drawGrid();
		input.hide()
		button.hide()
		slider.show()
		for(let i=0; i<slider.value() ;i++){
			let temp = factorPair(i);
			temp.forEach(e =>{
				plot(e)
				plot(e.reverse())
			})
		}
	
	}
	else if(radio.value() == 2){
		primeGrid()
		slider.hide()
		input.show()
		button.show()
		
		
		if(seq != []){
			for(e in seq){
				primePlot(seq[e][0], seq[e][1], seq[e][2], e)
				console.log(seq[e][2])
				
			}
		}
		
		
	}
	
}

function drawGrid(){
	inter = Math.floor((SQ-OFF)/ SIZE)
	lineLen = inter * (SIZE-1)
	for(let i = 0; i < SIZE; i++){
		line(OFF, (SQ-OFF) - (i * inter), OFF + lineLen, (SQ-OFF) - (i * inter));
		textSize(12)
		textAlign(RIGHT)
		text(i, OFF -10, (SQ-OFF) - (i * inter) + 3) //Dont mind the magic numbres
	}
	
	var top = (SQ-OFF) - ((SIZE-1) * inter)
	
	for(let i = 0; i < SIZE; i++){
		line(OFF + (i * inter), top, OFF + (i * inter), SQ-OFF);
		textSize(12)
		textAlign(CENTER, BOTTOM)
		text(i, OFF + (i * inter), SQ-OFF+22) //Dont mind the magic numbres
		
	}
}

function primeGrid(){
	
	inter = Math.floor((SQ-OFF)/ SIZE)
	lineLen = inter * (SIZE-1)
	i = 0
	fill(color(0,0,0))
	for(p in primes){
		
		line(OFF, (SQ-OFF) - (i * inter), OFF + lineLen, (SQ-OFF) - (i * inter));
		textSize(12)
		textAlign(RIGHT)
		text(primes[p], OFF -10, (SQ-OFF) - (i * inter) + 3) //Dont mind the magic numbres
		i++
	}
	
	var top = (SQ-OFF) - ((SIZE-1) * inter)
	i = 0
	for(p in primes){

		line(OFF + (i * inter), top, OFF + (i * inter), SQ-OFF);
		textSize(12)
		textAlign(CENTER, BOTTOM)
		text(primes[p], OFF + (i * inter), SQ-OFF+22) //Dont mind the magic numbres
		i++
	}
}

function primePlot(x, y, even, num){
	if(x<SIZE && y < SIZE){
		if(even)fill(color(0, 0, 255))
		else fill(color(255, 0, 0))
		circle(OFF + (x * inter), (SQ-OFF) - (y * inter), 10)
		text(num, OFF + (x * inter)-5, (SQ-OFF) - (y * inter)-5)
	}
}

function plot(x,y){
	if(x<SIZE && y < SIZE){
		circle(OFF + (x * inter), (SQ-OFF) - (y * inter), 10)
	}
	
}

function plot(list){
	if(list[0]<SIZE && list[1] < SIZE){
		circle(OFF + (list[0] * inter), (SQ-OFF) - (list[1] * inter), 10)
	}
	
}

function runSeq(){
	seq = []
	n = input.value()
	pair = distinctPrimes(input.value())
	while(n > 8){
		if(n %2 ==0) seq.push([primes.indexOf(pair[0]), primes.indexOf(pair[1]), 0])
		else seq.push([primes.indexOf(pair[0]), primes.indexOf(pair[1]), 1])
		n =sumDiv(n)
		pair = distinctPrimes(n)
	}

	
}
