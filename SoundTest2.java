import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;

import javax.imageio.ImageIO;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;



public class SoundTest2 {
	
	private static class Peak{

		public int frame;
		public int start;
		public int end;
		public int[] peak;
		public double value;
		
		public Peak(int frame, int start, int[] peak, int end, double value){
			this.frame = frame;
			this.start = start;
			this.end = end;
			this.peak = peak;
			this.value = value;
		}
	}
	
	public static void pcmToWav(int sampleRate, int channel, int bitsPerSample, byte[] pcmBytes, String filePath) throws IOException{
		FileOutputStream output = new FileOutputStream(filePath);
		output.write('R');
		output.write('I');
		output.write('F');
		output.write('F');
		int dataSize = 36 + pcmBytes.length;
		output.write(dataSize & 0xFF);
		output.write(dataSize >> 8 & 0xFF);
		output.write(dataSize >> 16 & 0xFF);
		output.write(dataSize >> 24 & 0xFF);
		output.write('W');
		output.write('A');
		output.write('V');
		output.write('E');
		output.write('f');
		output.write('m');
		output.write('t');
		output.write(' ');
		output.write(16);
		output.write(0);
		output.write(0);
		output.write(0);
		output.write(1);
		output.write(0);
		output.write(channel);
		output.write(0);
		output.write(sampleRate & 0xFF);
		output.write(sampleRate >> 8 & 0xFF);
		output.write(sampleRate >> 16 & 0xFF);
		output.write(sampleRate >> 24 & 0xFF);
		int byteRate = sampleRate * channel * bitsPerSample / 8;
		output.write(byteRate & 0xFF);
		output.write(byteRate >> 8 & 0xFF);
		output.write(byteRate >> 16 & 0xFF);
		output.write(byteRate >> 24 & 0xFF);
		int blockAlign = channel * bitsPerSample / 8;
		output.write(blockAlign);
		output.write(0);
		output.write(bitsPerSample);
		output.write(0);
		output.write('d');
		output.write('a');
		output.write('t');
		output.write('a');
		output.write(pcmBytes.length & 0xFF);
		output.write(pcmBytes.length >> 8 & 0xFF);
		output.write(pcmBytes.length >> 16 & 0xFF);
		output.write(pcmBytes.length >> 24 & 0xFF);
		output.write(pcmBytes);
		output.flush();
		output.close();
	}
	
	public static double equal_loudness(double frequency) {
		double pivot = 1000.;
    	double h_pivot = ((1037918.48 - pivot * pivot) * (1037918.48 - pivot * pivot) + 1080768.16 * pivot * pivot) / ((9837328 - pivot * pivot) * (9837328 - pivot * pivot) + 11723776 * pivot * pivot);
    	double n_pivot = (pivot / (6.8966888496476 * Math.pow(10, -5))) * Math.sqrt(h_pivot / ((pivot * pivot + 79919.29) * (pivot * pivot + 1345600)));
    	double h_freq = ((1037918.48 - frequency * frequency) * (1037918.48 - frequency * frequency) + 1080768.16 * frequency * frequency) / ((9837328 - frequency * frequency) * (9837328 - frequency * frequency) + 11723776 * frequency * frequency);
    	double n_freq = (frequency / (6.8966888496476 * Math.pow(10, -5))) * Math.sqrt(h_freq / ((frequency * frequency + 79919.29) * (frequency * frequency + 1345600)));
    	return Math.abs(n_pivot / n_freq);
	}
	
	public static class Mel {
		double frequency;
		double magnitude;
		double angle;
		
		public Mel(double frequency, double magnitude, double angle) {
			this.frequency = frequency;
			this.magnitude = magnitude;
			this.angle = angle;
		}
	}
	
	public static byte[] mel_to_pcm(int sampleRate, double frame_sec, List<List<Mel>> list) {
		//byte[] pcmBytesArray = new byte[frame_sec * list.size() * sampleRate * 2];
		
		int frame_samples = (int)(frame_sec * sampleRate);
		
		byte[] pcmBytesArray = new byte[frame_samples * list.size() * 2];
		
		for(int i = 0; i < list.size(); i++) {
			List<Mel> mel_list = list.get(i);
			double[] amp = new double[mel_list.size()];
			for(int j = 0; j < mel_list.size(); j++) {
				//amp[j] = equal_loudness(mel_list.get(j).frequency) / equal_loudness(30);
				amp[j] = 1;
			}
			double frame_max = 0;
			double[] sample = new double[frame_samples];
			for(int j = 0; j < frame_samples; j++) {
				double common = 2 * Math.PI * ((j + 1.) / sampleRate);
				for(int k = 0; k < mel_list.size(); k++) sample[j] += Math.sin(common * mel_list.get(k).frequency + mel_list.get(k).angle) * amp[k] * mel_list.get(k).magnitude;
				if(frame_max < Math.abs(sample[j])) frame_max = Math.abs(sample[j]);
			}
			
			for(int j = 0; j < frame_samples; j++) {
				pcmBytesArray[(i * frame_samples + j) * 2] = (byte)((int)(sample[j] / frame_max * 32767) & 0xFF);
            	pcmBytesArray[(i * frame_samples + j) * 2 + 1] = (byte)((int)(sample[j] / frame_max * 32767) >> 8 & 0xFF);
			}
			
			/*
			for(int j = 0; j < frame_samples; j++) {
				double common = 2 * Math.PI * ((j + 1.) / sampleRate);
				double sample = 0;
				double weight = 0;
				for(int k = 0; k < mel_list.size(); k++) {
					sample = (sample * weight + Math.sin(common * mel_list.get(k).frequency + mel_list.get(k).angle) * 32767 * amp[k] * mel_list.get(k).magnitude) / (weight + mel_list.get(k).magnitude);
					weight = Math.max(weight, mel_list.get(k).magnitude);
				}
				sample *= weight;
				pcmBytesArray[(i * frame_samples + j) * 2] = (byte)((int)sample & 0xFF);
            	pcmBytesArray[(i * frame_samples + j) * 2 + 1] = (byte)((int)sample >> 8 & 0xFF);
			}
			*/
			
			//System.out.println(i);
		}
		
		return pcmBytesArray;
	}
	
	static final double A0 = 27.5 / 2;

	/*
	public static class Prewitt {
		int r;
		int size;
		ArrayList<int[]> x_h;
		ArrayList<int[]> y_h;
		int[] x_h_moving_sum;
		int[] y_h_moving_sum;
		ArrayList<int[][]> rslt_list;
		ArrayList<double[][]> gradient_tmp;

		public Prewitt (int size){
			//size must be greater than 3
			this.size = size;
			r = size / 2;
		}
		public void Add(int frame, int[] grayscale, boolean last) {
			if(frame == 0) {
				x_h = new ArrayList<>();
				y_h = new ArrayList<>();
				x_h_moving_sum = new int[grayscale.length];
				y_h_moving_sum = new int[grayscale.length];
				rslt_list = new ArrayList<>();
				gradient_tmp = new ArrayList<>();
			}
			int[] x_v = new int[grayscale.length];
			int[] y_v = new int[grayscale.length];
			int x_v_moving_sum = 0;
			int y_v_moving_sum = 0;
			for(int i = -r; i <= r; i++) {
				x_v_moving_sum += grayscale[Math.max(0, i - 1)] * ((i > 0 ? 1 : 0) - (i < 0 ? 1 : 0));
				y_v_moving_sum += grayscale[Math.max(0, i - 1)];
			}
			for(int i = 0; i < grayscale.length; i++) {
				x_v_moving_sum += grayscale[Math.max(0, i - (r + 1))] - grayscale[Math.max(0, i - 1)] - grayscale[i] + grayscale[Math.min(grayscale.length - 1, i + r)];
				y_v_moving_sum += grayscale[Math.min(grayscale.length - 1, i + r)] - grayscale[Math.max(0, i - (r + 1))];
				x_v[i] = x_v_moving_sum;
				y_v[i] = y_v_moving_sum;
			}
			for(int i = frame + (frame == 0 ? -r : 0); i <= frame + (last ? r : 0); i++) {
				x_h.add(x_v);
				y_h.add(y_v);
				for(int j = 0; j < grayscale.length; j++) {
					x_h_moving_sum[j] += x_v[j];
					y_h_moving_sum[j] += y_v[j] * ((i > 0 ? 1 : 0) - (i < 0 ? 1 : 0));
				}
				if(i < r) continue;
				double[][] gradient = new double[grayscale.length][2];
				for(int j = 0; j < grayscale.length; j++) {
					if(x_h_moving_sum[j] == y_h_moving_sum[j]) gradient[j][0] = gradient[j][1] = 0;
					else {
						double rad = Math.atan2(y_h_moving_sum[j], x_h_moving_sum[j]);
		        		if(rad < 0) rad += Math.PI * 2;
		        		gradient[j][0] = rad;
		        		gradient[j][1] = Math.sqrt(x_h_moving_sum[j] * x_h_moving_sum[j] + y_h_moving_sum[j] * y_h_moving_sum[j]);
					}
				}
				
				
				gradient_tmp.add(gradient);
				if(gradient_tmp.size() == 3) {
					int[][] canny = new int[grayscale.length][2];
					for(int j = 1; j < grayscale.length - 1; j++) {
						
						/
						int sector = (int)Math.floor((gradient_tmp.get(1)[j][0] + Math.PI / 4) / (Math.PI / 2)) % 4;
						if((sector % 2 == 0 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(1)[j - 1][1], gradient_tmp.get(1)[j + 1][1])) ||
							(sector % 2 == 1 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j][1], gradient_tmp.get(2)[j][1]))
	        			) {
	        			/

						/
						int sector = (int)Math.floor((gradient_tmp.get(1)[j][0] + Math.PI / 8) / (Math.PI / 4)) % 8;
						if((sector % 4 == 0 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(1)[j - 1][1], gradient_tmp.get(1)[j + 1][1])) ||
							(sector % 4 == 1 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j - 1][1], gradient_tmp.get(2)[j + 1][1])) ||
							(sector % 4 == 2 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j][1], gradient_tmp.get(2)[j][1])) ||
							(sector % 4 == 3 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j + 1][1], gradient_tmp.get(2)[j - 1][1]))
	        			) {
	        			/

							//canny[j][0] = (int)Math.floor((gradient_tmp.get(1)[j][0] + Math.PI / 8) / (Math.PI / 4)) % 8;
							canny[j][0] = (int)Math.floor((gradient_tmp.get(1)[j][0] + Math.PI / 4) / (Math.PI / 2)) % 4;
							canny[j][1] = Math.min(0xFF, (int)gradient_tmp.get(1)[j][1]);
						//}
					}
					rslt_list.add(canny);
					gradient_tmp.remove(0);
				}

				int[] x_v1 = x_h.remove(0);
				int[] y_v1 = y_h.get(y_h.size() / 2);
				int[] y_v2 = y_h.get(y_h.size() / 2 + 1);
				int[] y_v3 = y_h.remove(0);
				for(int j = 0; j < grayscale.length; j++) {
					x_h_moving_sum[j] -= x_v1[j];
					y_h_moving_sum[j] += y_v3[j] - y_v1[j] - y_v2[j];
				}
			}
		}
	}
	*/
	
	
	
	public static void main(String[] args) {
		
		//Math.sin(3x + Math.PI / 4)
		//Math.sin(3x + 0)
		
		//System.out.println(Math.cos((1 - 1.1)) * Math.sin((1 + 1.1)));
		
		/*
		for(int i = 1; i < 1000; i++) {
			double freq1 = 1;
			double freq2 = freq1 * (i / 10.);
			double max = 0;
			for(int j = 0; j < 44100; j++) {
				double common = 2 * Math.PI * ((j + 1.) / 44100);
				double s = Math.sin(common * freq1) + Math.sin(common * freq2);
				if(max < Math.abs(s)) max = Math.abs(s);
			}
			//System.out.println(i + " : " + max);
		}
		*/
		
		double amp1 = 1;
		double amp2 = 1;
		double angle1 = Math.PI / 15;
		double angle2 = 0;
		
		//System.out.println(Math.sqrt(Math.pow(Math.cos(angle1) + Math.cos(angle2), 2) + Math.pow(Math.sin(angle1) + Math.sin(angle2), 2)));
		//System.out.println(Math.sqrt(amp1 + amp2 + 2 * amp1 * amp2 * Math.cos(angle1 - angle2)));
		
		//System.out.println(Math.sin(Math.PI));
		double freq1 = 3;
		double freq2 = 4;
		
		System.out.println(Math.sin(Math.PI / 2));
		
		/*
		System.out.println(Math.sin(3 * (Math.PI / 2)));
		System.out.println(Math.sin(3 * (Math.PI / 2 * 3)));
		
		System.out.println(3 * (Math.PI / 2) / Math.PI);
		
		System.out.println(Math.sin((Math.PI / 2)));
		System.out.println(Math.sin((Math.PI / 2 * 3)));
		*/
		
		try {
			long stamp = System.currentTimeMillis();
			int sampleRateOrigin = 44100;
			double sampleRate = sampleRateOrigin;

            double frame_sec = 0.02;
            //System.out.println(frame_sec);
            double win_sec = 0.1;
            int frame_len = (int)(sampleRate * frame_sec);
            double[] win = new double[(int)(sampleRate * win_sec)];
            double win_sum = 0;
            for(int i = 0; i < win.length; i++) {
            	double a0 = 0.21557895;
            	double a1 = 0.41663158;
            	double a2 = 0.277263158;
            	double a3 = 0.083578947;
            	double a4 = 0.006947368;
            	//win[i] = a0 - a1 * Math.cos(2.0 * Math.PI * i / (win.length - 1)) + a2 * Math.cos(4.0 * Math.PI * i / (win.length - 1)) - a3 * Math.cos(6.0 * Math.PI * i / (win.length - 1)) + a4 * Math.cos(8.0 * Math.PI * i / (win.length - 1));
            	
            	//win[i] = 0.54 - 0.46 * Math.cos(2.0 * Math.PI * i / (win.length - 1));
            	
            	//win[i] = 0.5 * (1 - Math.cos(2.0 * Math.PI * i / (win.length - 1)));
            	win[i] = 1;
            	
            	win_sum += win[i];
            }
            for(int i = 0; i < win.length; i++) win[i] *= 1 / win_sum;

            
            int n_fft = 2;
            int mel_compression = 2;
            
            while(n_fft < win.length) n_fft <<= 1;
			//while(n_fft < win.length || A0 * (Math.pow(2, (keyStart - 0.5 + 1. / mel_compression) / 12) - Math.pow(2, (keyStart - 0.5) / 12)) < sampleRate / n_fft) n_fft <<= 1;
            
            n_fft = 32768;
            
            System.out.println("1px : " + (sampleRate / n_fft) + "Hz");
            
            System.out.println("win : " + win.length + ", n_fft : " + n_fft);
            
            
            double[] sampleBuffer = new double[win.length];
            int sampleOffset = sampleBuffer.length - frame_len;
            
            double equal_loudness_max = 0;
            double[] eq = new double[n_fft / 2];
            for(int i = 0; i < eq.length; i++) {
            	double pivot = 1000.;
            	double freq = i;
            	double h_pivot = ((1037918.48 - pivot * pivot) * (1037918.48 - pivot * pivot) + 1080768.16 * pivot * pivot) / ((9837328 - pivot * pivot) * (9837328 - pivot * pivot) + 11723776 * pivot * pivot);
            	double n_pivot = (pivot / (6.8966888496476 * Math.pow(10, -5))) * Math.sqrt(h_pivot / ((pivot * pivot + 79919.29) * (pivot * pivot + 1345600)));
            	double h_freq = ((1037918.48 - freq * freq) * (1037918.48 - freq * freq) + 1080768.16 * freq * freq) / ((9837328 - freq * freq) * (9837328 - freq * freq) + 11723776 * freq * freq);
            	double n_freq = (freq / (6.8966888496476 * Math.pow(10, -5))) * Math.sqrt(h_freq / ((freq * freq + 79919.29) * (freq * freq + 1345600)));
            	eq[i] = Math.abs(n_freq / n_pivot);
            	if(equal_loudness_max < eq[i]) equal_loudness_max = eq[i];
            }

            for(int i = 0; i < eq.length; i++) eq[i] /= equal_loudness_max;
            BufferedImage img = new BufferedImage(eq.length, 500, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < img.getWidth(); i++) {
            	img.setRGB(i, (int)((img.getHeight() - 1) * (1 - eq[i])), 0xFFFFFFFF);
            }
            ImageIO.write(img, "PNG", new File("output_eq.png"));//
            
            

            int mel_comp = 20;
            int[][] mel_range = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp)][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key = (int)Math.round(Math.max(-1, 12 * Math.log((double)i * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	if(0 <= key && key < mel_range.length) {
            		if(mel_range[key][0] == 0) mel_range[key][0] = i;
            		mel_range[key][1] = i;
            	}
            }
            
            double bin_Hz = sampleRate / n_fft;
            
            /*
            int mel_comp = 10;
            int keyCount = 88;
            int keyStart = 0;
            int[][] mel_range = new int[keyCount * mel_comp][2];
            double[] mel_freq = new double[mel_range.length];
            double[][] mel_weight = new double[mel_range.length][];
            for(int i = 0; i < mel_range.length; i++) {
            	mel_freq[i] = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0;
            	double center = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0 / bin_Hz;
            	//double center = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0;
            	double left = center * Math.pow(2, -1. / 12);
            	double right = center * Math.pow(2, 1. / 12);
            	
            	//System.out.println(i + " : " + center + ", " + left + " ~ " + right);
            	
            	//if(i == 0) System.out.println(center - center * Math.pow(2, -0.1 / 12));
            	mel_range[i][0] = (int)Math.max(0, (int)Math.ceil(left));
				mel_range[i][1] = (int)Math.min(eq.length - 1, (int)Math.floor(right));
				mel_weight[i] = new double[mel_range[i][1] - mel_range[i][0] + 1];
				//double weight_sum = 0;
				for(int j = mel_range[i][0]; j <= mel_range[i][1]; j++) {
					double weight = 0;
					if(j < center) weight = (j - left) / (center - left);
					else weight = (j - right) / (center - right);
					mel_weight[i][j - mel_range[i][0]] = weight;

					//weight_sum += weight;
				}
				//for(int j = mel_range[i][0]; j <= mel_range[i][1]; j++) mel_weight[i][j - mel_range[i][0]] /= weight_sum;
            }
            */
            

            FileInputStream fis = new FileInputStream(new File("sample/drum.pcm"));
            
            ArrayList<double[]> mag_list = new ArrayList<>();
            ArrayList<double[]> phase_list = new ArrayList<>();
            ArrayList<double[]> freq_list = new ArrayList<>();
            
            int frame = 0;//임시
            
            int channelCount = 1;
            int div = channelCount * 2;
            byte[] buffer = new byte[4096];
            int bufferSize = 0;
            
            double mag_max = 0;

            
            ArrayList<ArrayList<Peak>> peaks_list = new ArrayList<>();
            
            FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            
            List<Double> frame_list = new ArrayList<>();
            
            List<List<Mel>> mel_frames = new ArrayList<>();

            
            boolean bb = false;
            
            while((bufferSize = fis.read(buffer, 0, buffer.length)) != -1) {
            	for(int i = 0; i < bufferSize / div; i++) {
            		int samp = 0;
            		for(int j = 0; j < channelCount; j++) samp += (buffer[i * div + j * 2 + 1] << 8) | (buffer[i * div + j * 2] & 0xFF);
        			sampleBuffer[sampleOffset++] = samp / channelCount / 32767.;
                    if(sampleOffset == sampleBuffer.length) {
                    	
                    	
                    	
                    	Complex[] complex = new Complex[n_fft];
                    	
                    	double rms = 0;
                    	
    					for(int j = 0; j < n_fft; j++) {
    						complex[j] = new Complex(j < sampleBuffer.length ? sampleBuffer[j] * Math.sqrt(win[j]) : 0, 0);
    						
    						if(j < sampleBuffer.length) rms += sampleBuffer[j] * sampleBuffer[j] * win[j];
    					}
    					
    					rms = Math.sqrt(rms / sampleBuffer.length);
    					
    					complex = fft.transform(complex, TransformType.FORWARD);

    					/*
    					double[] mel = new double[mel_range.length];
    					double[] freq = new double[mel_range.length];
    					double[] phase = new double[mel_range.length];
    					for(int j = 0; j < mel_range.length; j++) {
    						for(int k = mel_range[j][0]; k <= mel_range[j][1]; k++) {
    							if(rms < complex[k].abs() && mel[j] < complex[k].abs()) {
    								freq[j] = k * bin_Hz;
    								mel[j] = complex[k].abs();
    								phase[j] = Math.atan2(complex[k].getImaginary(), complex[k].getReal());
    							}
    						}
    						if(phase[j] < 0) phase[j] += 2 * Math.PI;
    						
    						
    						if(mag_max < mel[j]) mag_max = mel[j];
    					}
    					mag_list.add(mel);
    					phase_list.add(phase);
    					freq_list.add(freq);
    					*/

    					
    					
    					/*
    					double[] mel = new double[mel_range.length];
    					double[] freq = new double[mel_range.length];
    					double[] phase = new double[mel_range.length];
    					for(int j = 0; j < mel_range.length; j++) {
    						for(int k = mel_range[j][0]; k <= mel_range[j][1]; k++) {
    							mel[j] += complex[k].abs() * mel_weight[j][k - mel_range[j][0]];
    							freq[j] = k * bin_Hz;
    							phase[j] = Math.atan2(complex[k].getImaginary(), complex[k].getReal());
    						}
    						if(mag_max < mel[j]) mag_max = mel[j];
    					}
    					mag_list.add(mel);
    					phase_list.add(phase);
    					freq_list.add(freq);
    					*/

    					double[] mag = new double[n_fft / 2];
    					double max = 0;
    					for(int j = 0; j < mag.length; j++) {
    						mag[j] = complex[j].abs();
    						if(max < mag[j]) max = mag[j];
    					}

    					double[] value = mag;
    					ArrayList<Peak> peaks = new ArrayList<>();
						for(int j = 0; j < value.length; j++) {
							if(value[j] > 0) {
								int start = j;
								int[] peak = new int[2];
								peak[0] = j;
								while(j + 1 < value.length && value[j] <= value[j + 1]) {
									if(value[j] < value[j + 1]) peak[0] = j + 1;
									j++;
								}
								peak[1] = j;
								while(j + 1 < value.length && 0 < value[j + 1] && value[j + 1] <= value[j]) {
									j++;
								}
								int end = j;
								
								//double value = dbfs[peak[0]];

								//double value = Math.max(0, 20 * Math.log10(scal[peak[0]]) + dbNorm) / dbNorm;
								//if(valid[peak[0]])
								//if(end - start + 1 > mel_px)
								if(value[peak[0]] > max * 0.1) peaks.add(new Peak(frame, start, peak, end, value[peak[0]]));
								//sum += value;
								//if(comp_max < value) comp_max = value;
							}
						}
						peaks_list.add(peaks);
    					
    					
						/*
    					double[] mag = new double[n_fft / 2];
    					double[] freq = new double[n_fft / 2];
    					double[] phase = new double[n_fft / 2];
    					for(int j = 0; j < mag.length; j++) {
    						freq[j] = j * bin_Hz;
    						mag[j] = complex[j].abs();
    						phase[j] = Math.atan2(complex[j].getImaginary(), complex[j].getReal());
    						if(phase[j] < 0) phase[j] += 2 * Math.PI;
    					}
    					mag_list.add(mag);
    					freq_list.add(freq);
    					phase_list.add(phase);
    					*/


    					sampleOffset = sampleBuffer.length - frame_len;
                        System.arraycopy(sampleBuffer, frame_len, sampleBuffer, 0, sampleOffset);
                        frame++;
                        
                        if(frame == 1000) bb = true;
                    }
            	}
            	if(bb) break;
            }
            System.out.println(System.currentTimeMillis() - stamp);
            
            BufferedImage peak = new BufferedImage(peaks_list.size(), n_fft / 2, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < peaks_list.size(); i++) {
            	ArrayList<Peak> peaks = peaks_list.get(i);
            	for(int j = 0; j < peaks.size(); j++) {
            		for(int k = peaks.get(j).start; k <= peaks.get(j).end; k++) {
            			int col = 0;
            			
            			//int degree = (int)(peaks.get(j).value / (scale - 1.) * 0xFF);
            			int degree = (int)(peaks.get(j).value * 0xFF);
            			
            			if(peaks.get(j).peak[0] <= k && k <= peaks.get(j).peak[1]) col = 0xFFFF0000;
            			else if(k < peaks.get(j).peak[0]) col = (degree << 24 | degree << 8) ; //0xFF00FF00;
            			else if(k > peaks.get(j).peak[1]) col = (degree << 24 | degree);
            			peak.setRGB(i, n_fft / 2 - 1 - k, col);
            		}
            	}
            }
            ImageIO.write(peak, "PNG", new File("output_peak.png"));//
            
            BufferedImage dbfs = new BufferedImage(mag_list.size(), mag_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mag_list.size(); i++) {
            	double[] mel = mag_list.get(i);
            	for(int j = 1; j < mel.length - 1; j++) {
            		int col = (int)(mel[j] / mag_max * 0xFF);
            		dbfs.setRGB(i, mag_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mel[j] > Math.max(mel[j - 1], mel[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            
            
            
            int frame_samples = (int)(frame_sec * sampleRate);
            
            System.out.println("frame_size : " + mag_list.size());
            
            double[] samplesArray = new double[frame_samples * mag_list.size()];
            //byte[] pcmBytesArray = new byte[frame_samples * mag_list.size() * 2];
            double sampleMax = 0;
            for(int i = 0; i < mag_list.size(); i++) {
            	double[] mag = mag_list.get(i);
            	double[] phase = phase_list.get(i);
            	double[] freq = freq_list.get(i);
            	

            	for(int j = 0; j < frame_samples; j++) {
            		double common = 2 * Math.PI * (j + 0.) / sampleRate;
            		//for(int k = 1; k < mag.length - 1; k++) if(mag[k] > Math.max(mag[k - 1], mag[k + 1])) samplesArray[i * frame_samples + j] += Math.sin(common * freq[k] + phase[k]) * mag[k];
            		
            		for(int k = 0; k < mag.length; k++) samplesArray[i * frame_samples + j] += Math.sin(common * freq[k] + phase[k]) * mag[k];
            		
            		sampleMax = Math.max(sampleMax, Math.abs(samplesArray[i * frame_samples + j]));
            	}
            	
            	System.out.println(i);
            	
            	//for(int j = 0; k < mag.length; k++)
            	
            	/*
        		for(int i = 0; i < list.size(); i++) {
        			List<Mel> mel_list = list.get(i);
        			double[] amp = new double[mel_list.size()];
        			for(int j = 0; j < mel_list.size(); j++) {
        				//amp[j] = equal_loudness(mel_list.get(j).frequency) / equal_loudness(30);
        				amp[j] = 1;
        			}
        			double frame_max = 0;
        			double[] sample = new double[frame_samples];
        			for(int j = 0; j < frame_samples; j++) {
        				double common = 2 * Math.PI * ((j + 1.) / sampleRate);
        				for(int k = 0; k < mel_list.size(); k++) sample[j] += Math.sin(common * mel_list.get(k).frequency + mel_list.get(k).angle) * amp[k] * mel_list.get(k).magnitude;
        				if(frame_max < Math.abs(sample[j])) frame_max = Math.abs(sample[j]);
        			}
        			for(int j = 0; j < frame_samples; j++) {
        				pcmBytesArray[(i * frame_samples + j) * 2] = (byte)((int)(sample[j] / frame_max * 32767) & 0xFF);
                    	pcmBytesArray[(i * frame_samples + j) * 2 + 1] = (byte)((int)(sample[j] / frame_max * 32767) >> 8 & 0xFF);
        			}
        		}
        		*/
            }
            
            System.out.println(sampleMax);
            byte[] pcmBytesArray = new byte[samplesArray.length * 2];
            for(int i = 0; i < samplesArray.length; i++) {
            	pcmBytesArray[i * 2] = (byte)((int)(samplesArray[i] / sampleMax * 32767) & 0xFF);
            	pcmBytesArray[i * 2 + 1] = (byte)((int)(samplesArray[i] / sampleMax * 32767) >> 8 & 0xFF);
            }
            
            /*
            BufferedImage dbfs = new BufferedImage(dbfs_list.size(), dbfs_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < dbfs_list.size(); i++) {
            	double[] mel = dbfs_list.get(i);
            	double[] phase = dbfs2_list.get(i);
            	int[] freq = freq_list.get(i);
            	
            	List<Mel> mels = new ArrayList<>();
            	
            	for(int j = 1; j < mel.length - 1; j++) {
            		int col = (int)(mel[j] / mag_max * 0xFF);
            		dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		if(mel[j] > Math.max(mel[j - 1], mel[j + 1])) mels.add(new Mel(freq[j] * bin_Hz, mel[j] / mag_max, phase[j]));
            	}
            	mel_frames.add(mels);
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            */
            
            
            //byte[] pcmBytesArray = mel_to_pcm((int)sampleRate, frame_sec, mel_frames);
            
            pcmToWav((int)sampleRate, channelCount, 16, pcmBytesArray, "output.wav");

            
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
}
