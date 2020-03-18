// Start up server
s.boot;

// Read in the sound file
//b = Buffer.read(s, "C:/Users/Ilhan/Desktop/rhythm_2.wav");
b = Buffer.read(s, "C:/Users/Ilhan/Desktop/woef.wav");

// Make a buffer to store feature vector values (energy method via Sethares)
f = Array.new(50);

// Sample some overlapping segments and store the total energy
for (0,49, {
	arg i;
	f.add(0.0);
	for (i*10, i*10 + 100, {
		arg j;
		b.get(j, {
			arg value;
			f.put(i, f[i] + value);
		});
	});
});

// Show energy method array contents
f.postln;

// Plot energy array as graph
f.plot;

// Plot the music
b.plot;

// Play the music
b.play;

// FFT plotting //

c = Buffer.alloc(s,2048,1);

// Create a playback of the sound file
(
x = { var in, chain, chainB, chainC;
	in = PlayBuf.ar(1,b.bufnum);
	//in = WhiteNoise.ar;
    chain = FFT(c, in);
	Pan2.ar(IFFT(chain));
    //0.01 * Pan2.ar(IFFT(chain));
}.play(s);
)

// Plot the magnitude of the sound progressively
(
Routine({
    2.do{arg i;
        c.getToFloatArray(action: { arg array;
            var z, x;
            z = array.clump(2).flop;
            // Initially data is in complex form
            z = [Signal.newFrom(z[0]), Signal.newFrom(z[1])];
            x = Complex(z[0], z[1]);

			p = x.phase;
			q = p;
			p = p.do({
				// "Unwrap" each phase element to preserve continuity
				arg item, j;
				if ((item < 3.14159) && (item > 0), { q.put(j, item); });
				if (item < 0, { q.put(j, (item + 3.14159)) });
				if (item > 3.14159, { q.put(j, (item - 3.14159)) });
			});

            { q.plot('Initial', Rect(200, 600-(200*i), 700, 200)) }.defer
        });
        0.0001.wait;
}}).play
)

// Clear the playback object
x.free;

// Clear buffers
(
b.free;
f.free;
t.free;
)

// High frequency sinusoid for possible beat tracking marker
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

k = Pbind(\instrument, \pn, \delta, 1, \preamp, 4.5, \amp, 0.32).play(quant: 1);

// Stop the sounds
p.stop;
k.stop;

// Shut down the server
s.quit;