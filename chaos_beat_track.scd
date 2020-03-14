s.boot;
// reads in the sound file
b = Buffer.read(s, "C:/Users/Ilhan/Desktop/rhythm_2.wav");

// Make a buffer to store feature vector values
f = Array.new(100);
//g = Array.new(5000);

for (0,9, {
	arg i;
	f.add(0.0);
	for (i*2, i*2 + 10, {
		arg j;
		//if ((j == (i*2)), { f.put(i,0.0) });
		b.get(j, {
			arg value;
			f.put(i, f[i] + value);
		});
	});
});

f.postln;
g.postln;

f.plot;

b.plot;
// play the contents
b.play;
// Clear the buffer
b.free;
f.free;

// totally cheesy, but who could resist?
(
SynthDef(\kik, { |out, preamp = 1, amp = 1|
    var //   freq = EnvGen.kr(Env([400, 66], [0.08], -3)),
        sig = SinOsc.ar(10000, 0.5pi, preamp).distort * amp
            * EnvGen.kr(Env([0, 1, 0.8, 0], [0.01, 0.1, 0.2]), doneAction: Done.freeSelf);
    Out.ar(out, sig ! 2);
}).add;

(
SynthDef(\pn, { |out, preamp = 1, amp = 1|
    var //   freq = EnvGen.kr(Env([400, 66], [0.08], -3)),
       /* sig = PinkNoise.ar(0.2).distort * amp
            * EnvGen.kr(Env([0, 1, 0.8, 0], [0.01, 0.1, 0.2]), doneAction: Done.freeSelf);*/
	sig = SinOsc.ar(5000,0,0.2) * amp
	* EnvGen.kr(Env([0, 1, 0.8, 0], [0.01, 0.1, 0.2]), doneAction: Done.freeSelf);
    Out.ar(out, sig ! 2);
}).add;

// before you play:
// what do you anticipate '\delta, 1' will do?
k = Pbind(\instrument, \pn, \delta, 1, \preamp, 4.5, \amp, 0.32).play(quant: 1);

p.stop;
k.stop;

s.quit;
