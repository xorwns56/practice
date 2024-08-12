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



public class SoundTest {
	private static class IndexValuePair{
		public int index;
		public double value;
		public IndexValuePair(int index, double value) {
			this.index = index;
			this.value = value;
		}
	}
	
	private static class Peak{
		public int frame;
		public int start;
		public int end;
		public int[] peak;
		public double value;
		public double phase;
		public double freq;
		
		public Peak(int frame, int start, int[] peak, int end, double value, double phase, double freq){
			this.frame = frame;
			this.start = start;
			this.end = end;
			this.peak = peak;
			this.value = value;
			this.phase = phase;
			this.freq = freq;
		}
	}

	
	private static class Data{
		public double freq;
		public double phase;
		public double value;
	}
	
	public static class Sobel {
		ArrayList<double[]> x_h;
		ArrayList<double[]> y_h;
		ArrayList<double[][]> rslt_list;
		ArrayList<double[][]> gradient_tmp;
		public void Add(int frame, double[] value, boolean last) {
			if(frame == 0) {
				x_h = new ArrayList<>();
				y_h = new ArrayList<>();
				rslt_list = new ArrayList<>();
				gradient_tmp = new ArrayList<>();
			}
			double[] x_v = new double[value.length];
			double[] y_v = new double[value.length];
			for(int i = 0; i < value.length; i++) {
				x_v[i] = value[Math.min(value.length - 1, i + 1)] - value[Math.max(0, i - 1)];
				y_v[i] = value[Math.max(0, i - 1)] + value[i] * 2 + value[Math.min(value.length - 1, i + 1)];
			}
			for(int i = frame + (frame == 0 ? -1 : 0); i <= frame + (last ? 1 : 0); i++) {
				x_h.add(x_v);
				y_h.add(y_v);
				if(i < 1) continue;
				double[][] gradient = new double[value.length][2];
				for(int j = 0; j < value.length; j++) {
					double x_sum = x_h.get(0)[j] + x_h.get(1)[j] * 2 + x_h.get(2)[j];
					double y_sum = y_h.get(2)[j] - y_h.get(0)[j];
					if(x_sum == y_sum) gradient[j][0] = gradient[j][1] = 0;
					else {
						double rad = Math.atan2(y_sum, x_sum);
		        		if(rad < 0) rad += Math.PI * 2;
		        		
		        		
		        		gradient[j][0] = rad;
		        		gradient[j][1] = Math.sqrt(x_sum * x_sum + y_sum * y_sum);
					}
				}
				x_h.remove(0);
				y_h.remove(0);
				rslt_list.add(gradient);
				/*
				gradient_tmp.add(gradient);
				if(gradient_tmp.size() < 2) continue;
				double[][] canny = new double[value.length][2];
				if(gradient_tmp.size() == 3) {
					for(int j = 1; j < value.length - 1; j++) {
						
						/
						double a = 160. / 180. * Math.PI;
	        			double angle = gradient_tmp.get(1)[j][0] + a / 2;
	        			int sector = (int)(Math.floor(angle / Math.PI) * 2 + (Math.floor(angle % Math.PI / a) > 0 ? 1 : 0)) % 4;
						if((sector % 2 == 0 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(1)[j - 1][1], gradient_tmp.get(1)[j + 1][1])) ||
							(sector % 2 == 1 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j][1], gradient_tmp.get(2)[j][1]))
	        			) {
							canny[j][0] = gradient_tmp.get(1)[j][0];
							canny[j][1] = Math.min(1, gradient_tmp.get(1)[j][1]);
						}
						/

						int sector = (int)Math.floor((gradient_tmp.get(1)[j][0] + Math.PI / 8) / (Math.PI / 4)) % 8;
						if((sector % 4 == 0 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(1)[j - 1][1], gradient_tmp.get(1)[j + 1][1])) ||
							(sector % 4 == 1 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j - 1][1], gradient_tmp.get(2)[j + 1][1])) ||
							(sector % 4 == 2 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j][1], gradient_tmp.get(2)[j][1])) ||
							(sector % 4 == 3 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j + 1][1], gradient_tmp.get(2)[j - 1][1]))
	        			) {
							canny[j][0] = gradient_tmp.get(1)[j][0];
							canny[j][1] = Math.min(1, gradient_tmp.get(1)[j][1]);
						}


					}
					gradient_tmp.remove(0);
				}
				rslt_list.add(canny);
				*/

			}
		}
	}
	
	/*
	private static class Peak{
		public int[] freq_index;
		public int frame;
		public int start;
		public int end;
		public int[] peak;
		public double value;
		public double phase;
		public double freq;
		
		public Peak(int frame, int start, int[] peak, int end, double value, double phase, double freq){
			this.frame = frame;
			this.start = start;
			this.end = end;
			this.peak = peak;
			this.value = value;
			this.phase = phase;
			this.freq = freq;
		}
	}
	*/
	/*
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
	*/
	
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
	
	public static class Wave {
		double frequency;
		double amplitude;
		double angle;
		
		public Wave(double frequency, double amplitude, double angle) {
			this.frequency = frequency;
			this.amplitude = amplitude;
			this.angle = angle;
		}
	}
	
	public static byte[] mel_to_pcm(int sampleRate, double frame_sec, List<List<Wave>> list) {
		//byte[] pcmBytesArray = new byte[frame_sec * list.size() * sampleRate * 2];
		
		int frame_samples = (int)(frame_sec * sampleRate);
		
		byte[] pcmBytesArray = new byte[frame_samples * list.size() * 2];
		
		for(int i = 0; i < list.size(); i++) {
			List<Wave> mel_list = list.get(i);
			double[] amp = new double[mel_list.size()];
			for(int j = 0; j < mel_list.size(); j++) {
				//amp[j] = equal_loudness(mel_list.get(j).frequency) / equal_loudness(30);
				amp[j] = 1;
			}
			double frame_max = 0;
			double[] sample = new double[frame_samples];
			for(int j = 0; j < frame_samples; j++) {
				double common = 2 * Math.PI * ((j + 1.) / sampleRate);
				for(int k = 0; k < mel_list.size(); k++) sample[j] += Math.sin(common * mel_list.get(k).frequency + mel_list.get(k).angle) * amp[k] * mel_list.get(k).amplitude;
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
	
	static final double A0 = 27.5;

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
			FileInputStream fis = new FileInputStream(new File("sample/sampleTest.pcm"));
			int sampleRateOrigin = 44100;
			double sampleRate = sampleRateOrigin;

            double frame_sec = 0.02;
            //System.out.println(frame_sec);
            double win_sec = frame_sec * 2;
            //double win_sec = 0.02;
            int frame_len = (int)(sampleRate * frame_sec);
            double[] rect_win = new double[(int)(sampleRate * win_sec)];
            double[] hann_win = new double[rect_win.length];
            double rect_win_sum = 0;
            double hann_win_sum = 0;
            for(int i = 0; i < rect_win.length; i++) {
            	double a0 = 0.21557895;
            	double a1 = 0.41663158;
            	double a2 = 0.277263158;
            	double a3 = 0.083578947;
            	double a4 = 0.006947368;
            	//win[i] = a0 - a1 * Math.cos(2.0 * Math.PI * i / (win.length - 1)) + a2 * Math.cos(4.0 * Math.PI * i / (win.length - 1)) - a3 * Math.cos(6.0 * Math.PI * i / (win.length - 1)) + a4 * Math.cos(8.0 * Math.PI * i / (win.length - 1));
            	
            	//win[i] = 0.54 - 0.46 * Math.cos(2.0 * Math.PI * i / (win.length - 1));
            	rect_win[i] = 1;
            	hann_win[i] = 0.5 * (1 - Math.cos(2.0 * Math.PI * i / (hann_win.length - 1)));
            	
            	rect_win_sum += rect_win[i];
            	hann_win_sum += hann_win[i];
            }
            for(int i = 0; i < rect_win.length; i++) {
            	rect_win[i] *= 1 / rect_win_sum;
            	hann_win[i] *= 1 / hann_win_sum;
            }

            
            int n_fft = 2;
            
            while(n_fft < rect_win.length) n_fft <<= 1;
			//while(n_fft < win.length || A0 * (Math.pow(2, (keyStart - 0.5 + 1. / mel_compression) / 12) - Math.pow(2, (keyStart - 0.5) / 12)) < sampleRate / n_fft) n_fft <<= 1;

            //n_fft *= 2 * 2;
            n_fft = 32768;
            
            System.out.println("1px : " + (sampleRate / n_fft) + "Hz");
            
            System.out.println("win : " + rect_win.length + ", n_fft : " + n_fft);
            
            
            double[] sampleBuffer = new double[rect_win.length];
            int sampleOffset = sampleBuffer.length - frame_len;
            
            /*
            int mel_comp = 20;
            int[][] comp_range = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp) + 1][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key = (int)Math.round(Math.max(-1, 12 * Math.log((i + 0.5) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	
            	if(0 <= key) {
            		
            	}
            	
            	//System.out.println(i + " : " + key);
            }
            */

			

			int mel_comp = 20;
            int[][] comp_range = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp) + 1][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key1 = (int)Math.round(Math.max(-1, 12 * Math.log((i + 0.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	int key2 = (int)Math.round(Math.max(-1, 12 * Math.log((i + 1.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	for(int j = Math.max(0, key1); j <= key2; j++) {
            		if(comp_range[j][0] == 0) comp_range[j][0] = i;
            		comp_range[j][1] = i;
            	}
            }
            double[] comp_freq = new double[comp_range.length];
            for(int i = 0; i < comp_range.length; i++) {
            	//System.out.println(i + " : " + Math.pow(A0 ))
            	comp_freq[i] = Math.pow(2, ((double)i / mel_comp) / 12) * A0;
            }

            
            double bin_Hz = sampleRate / n_fft;
            

            double a_freq = 20;
            double b_freq = 20000;
            int a_f = (int)Math.ceil((12 * Math.log(a_freq / A0) / Math.log(2)) * mel_comp);
            int b_f = (int)Math.floor((12 * Math.log(Math.min(b_freq, sampleRate * 0.5) / A0) / Math.log(2)) * mel_comp);
            int[][] mel_range = new int[b_f - a_f + 1][2];
            double[] mel_freq = new double[mel_range.length];
            double[] mel_lerp = new double[mel_range.length];
            double[] mel_eq = new double[mel_range.length];
            double[][] mel_weight = new double[mel_range.length][];
            for(int i = 0; i < mel_range.length; i++) {
            	mel_freq[i] = Math.pow(2, ((i + a_f) / (double)mel_comp / 12)) * A0;
            	mel_lerp[i] = Math.pow(2, ((i + a_f) / (double)mel_comp / 12)) * A0 / bin_Hz;
            	//double center = mel_freq[i] / bin_Hz;
            	
            	//double center = (int)Math.round(Math.max(-1, 12 * Math.log((i + 0.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);

            	
            	
            	//double center = Math.pow(2, ((i - 0 - mel_comp / 2) / (double)mel_comp) / 12) * A0;
            	
            	
            	//System.out.println(center);
            	System.out.println(mel_lerp[i] * bin_Hz);
            	
            	
            	/*
            	double left = center * Math.pow(2, -influence / 12);
            	double right = center * Math.pow(2, influence / 12);
            	
            	//System.out.println(i + " : " + center + ", " + left + " ~ " + right);
            	
            	mel_range[i][0] = (int)Math.max(0, (int)Math.ceil(left));
				mel_range[i][1] = (int)Math.min(n_fft / 2, (int)Math.floor(right));
            	
				//System.out.println(i + " : " + center + ", " + mel_range[i][0] + " ~ " + mel_range[i][1]);

            	//if(i == 0) System.out.println(center - center * Math.pow(2, -0.1 / 12));
				mel_weight[i] = new double[mel_range[i][1] - mel_range[i][0] + 1];
				double weight_sum = 0;
				for(int j = mel_range[i][0]; j <= mel_range[i][1]; j++) {
					double weight = 0;
					if(j < center) weight = (j - left) / (center - left);
					else weight = (j - right) / (center - right);
					mel_weight[i][j - mel_range[i][0]] = weight;
					
					weight_sum += weight;
				}
				*/
            }
            
            
            int test_start = 13;
            int test_end = 13;
            
            List<Peak> test_peaks = new ArrayList<>();
            test_peaks.add(new Peak(0, 2, new int[] { 2, 2 }, 3, 0, 0, 0));
            test_peaks.add(new Peak(0, 4, new int[] { 5, 5 }, 6, 0, 0, 0));
            test_peaks.add(new Peak(0, 7, new int[] { 7, 7 }, 7, 0, 0, 0));
            test_peaks.add(new Peak(0, 8, new int[] { 9, 9 }, 10, 0, 0, 0));
            test_peaks.add(new Peak(0, 11, new int[] { 11, 11 }, 12, 0, 0, 0));
            int low = 0;
            int high = test_peaks.size();
            
            while(low < high) {
            	int mid = (low + high) / 2;
            	if(test_start <= test_peaks.get(mid).end) high = mid;
            	else low = mid + 1;
            }
            

            while(low < test_peaks.size() && test_peaks.get(low).start <= test_end) {
            	
            	System.out.println(test_peaks.get(low).start + ", " + test_peaks.get(low).end);
            	
            	low++;
            }
            
            
            double equal_loudness_max = 0;
            double[] eq = new double[n_fft / 2];
            for(int i = 0; i < eq.length; i++) {
				eq[i] = equal_loudness(i * bin_Hz);
            	if(equal_loudness_max < eq[i]) equal_loudness_max = eq[i];
            }
            

            
            PriorityQueue<IndexValuePair> pq = new PriorityQueue<>((pair1, pair2)->{ return Double.valueOf(pair2.value).compareTo(pair1.value); });
            int[][] mag_range = new int[eq.length][2];
            double mag_radius = 0.5;
            for(int j = 0; j < eq.length; j++) {
            	double fr1 = (j + 0.) * sampleRate / n_fft;
            	double fr2 = (j + 1.) * sampleRate / n_fft;
            	
            	//mag_range[j][0] = Math.max(0, j - mel_px / 2);
            	//mag_range[j][1] = Math.min(eq.length - 1, j + mel_px / 2);
				mag_range[j][0] = Math.max(0, (int)Math.floor(fr1 * Math.pow(2, -mag_radius / 12) * n_fft / sampleRate));
				mag_range[j][1] = Math.min(eq.length - 1, (int)Math.floor(fr2 * Math.pow(2, mag_radius / 12) * n_fft / sampleRate));
            }
            
            
            
            /*
            for(int i = 0; i < mel_eq.length; i++) mel_eq[i] /= equal_loudness_max;
            BufferedImage img = new BufferedImage(mel_eq.length, 500, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < img.getWidth(); i++) {
            	img.setRGB(i, (int)((img.getHeight() - 1) * (1 - mel_eq[i])), 0xFFFFFFFF);
            }
            ImageIO.write(img, "PNG", new File("output_eq.png"));//
            */

            ArrayList<double[]> sum_of_harmonic_list = new ArrayList<>();
            ArrayList<Integer> max_idx_of_harmonic_list = new ArrayList<>();
            double max_sum_of_harmonic = 0;
            int overtone_count = 10;
            double overtone_coef = 1.0;

            ArrayList<double[]> mel_list = new ArrayList<>();
            ArrayList<double[]> mag_eq_list = new ArrayList<>();
            ArrayList<double[]> mag_pq_list = new ArrayList<>();
            ArrayList<double[]> mag_list = new ArrayList<>();
            ArrayList<double[]> phase_list = new ArrayList<>();
            ArrayList<double[]> freq_list = new ArrayList<>();
            ArrayList<Complex[]> complex_list = new ArrayList<>();
            ArrayList<int[]> max_list = new ArrayList<>();
            
            Sobel sobel = new Sobel();
            ArrayList<double[][]> gradient_list = new ArrayList<>();
            
            int frame = 0;//임시
            
            int channelCount = 1;
            int div = channelCount * 2;
            byte[] buffer = new byte[4096];
            int bufferSize = 0;
            
            double mag_max = 0;
            double mag_eq_max = 0;
            double mag_pq_max = 0;
            
            double mel_max = 0;

            
            double gradient_max = 0;
            
            ArrayList<ArrayList<Peak>> peaks_list = new ArrayList<>();
            
            FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            
            List<Double> frame_list = new ArrayList<>();
            
            List<Wave[]> wave_list = new ArrayList<>();
            
            
            
            List<boolean[]> valid_list = new ArrayList<>();
            
            List<List<Wave>> mel_frames = new ArrayList<>();
            
            List<Double> peak_sum_list = new ArrayList<>();
            double max_peak_sum = 0;

            
            List<Peak> last_peaks = new ArrayList<>();
            //System.out.println(0.31198 / Math.PI);
            
            boolean bb = false;
            
            double[] prev_phase = new double[n_fft / 2];
            double[] prev_freq = new double[n_fft / 2];
            
            while((bufferSize = fis.read(buffer, 0, buffer.length)) != -1) {
            	for(int i = 0; i < bufferSize / div; i++) {
            		int samp = 0;
            		for(int j = 0; j < channelCount; j++) samp += (buffer[i * div + j * 2 + 1] << 8) | (buffer[i * div + j * 2] & 0xFF);
        			sampleBuffer[sampleOffset++] = samp / channelCount / 32767.;
                    if(sampleOffset == sampleBuffer.length) {

                    	Complex[] complex = new Complex[n_fft];
    					for(int j = 0; j < n_fft; j++) complex[j] = new Complex(j < sampleBuffer.length ? sampleBuffer[j] * rect_win[j] : 0, 0);
    					complex = fft.transform(complex, TransformType.FORWARD);
    					
    					//Wave[] wave = new Wave[n_fft / 2];
    					//for(int j = 0; j < wave.length; j++) wave[j] = new Wave(j * bin_Hz, complex[j].abs(), Math.atan2(complex[j].getImaginary(), complex[j].getReal()));
    					//wave_list.add(wave);
    					
    					double[] mag = new double[n_fft / 2];
    					double[] mag_eq = new double[mag.length];
    					
    					double[] phase = new double[mag.length];
    					
    					for(int j = 0; j < mag.length; j++) {
    						
    						phase[j] = Math.atan2(complex[j].getImaginary(), complex[j].getReal());
    						if(phase[j] < 0) phase[j] += 2 * Math.PI;
    						
    						
    						mag[j] = complex[j].abs();
    						if(mag_max < mag[j]) mag_max = mag[j];
    						mag_eq[j] = mag[j] / eq[j];
    						if(mag_eq_max < mag_eq[j]) mag_eq_max = mag_eq[j];
    					}
    					//mag_list.add(mag);
    					//mag_eq_list.add(mag_eq);
    					
    					phase_list.add(phase);
    					

    					double[] mag_pq = new double[mag.length];
    					int[] curr_range = new int[] {0, 0};
			            pq.add(new IndexValuePair(0, mag[0]));
						for(int j = 0; j < eq.length; j++) {
							//System.out.println(mag_range[j][0] + " ~ " + mag_range[j][1]);
							for(int k = curr_range[1] + 1; k <= mag_range[j][1]; k++) pq.add(new IndexValuePair(k, mag_eq[k]));
							curr_range = mag_range[j];
							while(pq.peek().index < curr_range[0]) pq.remove();
							mag_pq[j] = pq.peek().value;
							
							if(mag_pq_max < mag_pq[j]) mag_pq_max = mag_pq[j];
						}
						pq.clear();
    					
						mag_pq_list.add(mag_pq);

    					

    					double[] comp = new double[comp_range.length];
    					double[] comp2 = new double[comp_range.length];
    					for(int j = 0; j < comp.length; j++) {
    						for(int k = comp_range[j][0]; k <= comp_range[j][1]; k++) {
    							if(comp[j] < mag_eq[k]) {
    								comp[j] = mag_eq[k];
    								comp2[j] = mag[k];
    							}
    						}
    					}
    					//comp = mag_eq;
    					mag_list.add(comp2);
    					mag_eq_list.add(comp);

    					
    					/*
    					double[] mel = new double[mel_lerp.length];
						for(int j = 0; j < mel_lerp.length; j++) {
							int l = (int)mel_lerp[j];
							double mu = mel_lerp[j] - l;
							double mu2 = mu * mu;
							
							double y3 = l + 2 < mag.length ? mag_eq[l + 2] : 0;
							double y2 = l + 1 < mag.length ? mag_eq[l + 1] : 0;
							double y1 = mag_eq[l];
							double y0 = l - 1 < 0 ? 0 : mag_eq[l - 1];
							
							double a0 = y3 - y2 - y0 + y1;
							double a1 = y0 - y1 - a0;
							double a2 = y2 - y0;
							double a3 = y1;
							//dbfs[j] = Math.min(1, Math.max(0, a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3));
							//dbfs[j] = Math.min(1, Math.max(0, a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3));
							mel[j] = a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3;
							
							if(mel_max < mel[j]) mel_max = mel[j];
						}
						mel_list.add(mel);
						*/

    					
    					double peak_sum = 0;
    					
    					
    					for(int j = 0; j < mag_eq.length; j++) {
    						//peak_sum += mag_eq[j];
    						//peak_sum = Math.max(peak_sum, mag_eq[j]);
    					}
    					
    					
    					//comp = mag_eq;
    					
    					ArrayList<Peak> peaks = new ArrayList<>();
						for(int j = 0; j < comp.length; j++) {
							if(comp[j] > 0.00) {
								int start = j;
								int[] peak = new int[2];
								peak[0] = j;
								while(j + 1 < comp.length && comp[j] <= comp[j + 1]) {
									if(comp[j] < comp[j + 1]) peak[0] = j + 1;
									j++;
								}
								peak[1] = j;
								while(j + 1 < comp.length && 0.00 < comp[j + 1] && comp[j + 1] <= comp[j]) {
									j++;
								}
								int end = j;
								
								//lastPeak에서 연결되는 부분이 있는지 찾아야한다. 일단, start와 end에 해당 되는 모든 부분을 받아와야돼. start를 기준으로해서 이진검색 후 lastpeaks에서 반복돌린다.
								//그 후에 그 부분 들 중 소리가 비슷한 부분을 찾고(peak가 연결되면 확정시키나)?
								peaks.add(new Peak(frame, start, peak, end, comp[peak[0]], 0, 0));
								
								//peak_sum += comp[peak[0]];
								peak_sum = Math.max(peak_sum, comp[peak[0]]);
								
								//for(int k = comp_range[peak[0]][0]; k <= comp_range[peak[1]][1]; k++) valid[k] = true;
								
								//if(comp[peak[0]] > 0.001)
								//for(int k = comp_range[start][0]; k <= comp_range[end][1]; k++) valid[k] = true;
								//for(int k = max_idx[start]; k <= max_idx[end]; k++) valid[k] = true;
							}
						}
						peaks_list.add(peaks);
						last_peaks = peaks;
						//System.out.println(last_peaks.size());
						peak_sum_list.add(peak_sum);
						if(max_peak_sum < peak_sum) max_peak_sum = peak_sum;
						
						
						double[] sum_of_harmonic = new double[comp.length];
    					
    					int max_idx_of_harmonic = 0;
    					
						for(int k = 0; k < comp.length; k++) {
							for(int j = 0; j < overtone_count; j++) {
		                		int l = k + (int)Math.round(12 * mel_comp * Math.log(1 + j * overtone_coef) / Math.log(2));
		                		
								if(l < comp.length) sum_of_harmonic[k] += comp[l]; //comp[l];
		                    }
							if(sum_of_harmonic[max_idx_of_harmonic] < sum_of_harmonic[k]) max_idx_of_harmonic = k;
						}
						if(max_sum_of_harmonic < sum_of_harmonic[max_idx_of_harmonic]) max_sum_of_harmonic = sum_of_harmonic[max_idx_of_harmonic];
						max_idx_of_harmonic_list.add(max_idx_of_harmonic);
						sum_of_harmonic_list.add(sum_of_harmonic);

    					/*
    					double[] mel = new double[mel_range.length];
    					for(int j = 0; j < mel.length; j++) {
    						for(int k = mel_range[j][0]; k <= mel_range[j][1]; k++) {
    							mel[j] += mag[k] * mel_weight[j][k - mel_range[j][0]];
    							//comp[j] = Math.max(comp[j], wave[k].amplitude);
    						}
    					}
    					*/
    					
						/*
    					sobel.Add(frame, comp, false);
    					while(sobel.rslt_list.size() > 0) {
    						double[][] gradient = sobel.rslt_list.remove(0);
    						for(int j = 0; j < gradient.length; j++) {
    							if(gradient_max < gradient[j][1]) gradient_max = gradient[j][1];
    						}
    						gradient_list.add(gradient);
    					}
    					*/


    					sobel.Add(frame, mag_eq, false);
    					while(sobel.rslt_list.size() > 0) {
    						double[][] gradient = sobel.rslt_list.remove(0);
    						/*
    						for(int j = 0; j < gradient.length; j++) {
    							if(gradient_max < gradient[j][1]) gradient_max = gradient[j][1];
    						}
    						*/
    						
    						double[][] comp_gradient = new double[comp_range.length][2];
        					for(int j = 0; j < comp_gradient.length; j++) {
        						for(int k = comp_range[j][0]; k <= comp_range[j][1]; k++) {
        							if(comp_gradient[j][1] < gradient[k][1]) comp_gradient[j] = gradient[k];
        						}
        						
        						if(gradient_max < comp_gradient[j][1]) gradient_max = comp_gradient[j][1];
        					}
    						
    						
    						gradient_list.add(comp_gradient);
    					}
    					
    					/*
    					double[] mel = new double[comp_range.length];
    					int[] max_idx = new int[comp_range.length];
    					for(int j = 0; j < mel.length; j++) {
    						for(int k = comp_range[j][0]; k <= comp_range[j][1]; k++) {
    							
    							if(mel[j] < mag[k] / eq[k]) {
    								mel[j] = mag[k] / eq[k];
    								max_idx[j] = k;
    							}
    							
    							//mel[j] = Math.max(mel[j], wave[k].amplitude / eq[k]);
    						}
    						if(mag_max < mel[j]) mag_max = mel[j];
    					}
    					mag_list.add(mag);
    					mel_list.add(mel);
    					max_list.add(max_idx);
    					
						*/
    					
    					/*
    					boolean[] valid = new boolean[wave.length];
    					for(int j = 0; j < wave.length; j++) {
    						if(wave[j].amplitude * eq[j] > 0.01) valid[j] = true;
    					}
    					/
    					for(int j = 0; j < max_idx.length; j++) {
    						//if(wave[max_idx[j]].amplitude * eq[max_idx[j]] > 0.005)
    						valid[max_idx[j]] = true;
    					}
    					/
    					valid_list.add(valid);
    					*/


    					/*
    					ArrayList<Peak> peaks = new ArrayList<>();
						for(int j = 0; j < comp.length; j++) {
							if(comp[j] > 0.00) {
								int start = j;
								int[] peak = new int[2];
								peak[0] = j;
								while(j + 1 < comp.length && comp[j] <= comp[j + 1]) {
									if(comp[j] < comp[j + 1]) peak[0] = j + 1;
									j++;
								}
								peak[1] = j;
								while(j + 1 < comp.length && 0.00 < comp[j + 1] && comp[j + 1] <= comp[j]) {
									j++;
								}
								int end = j;
								
								peaks.add(new Peak(frame, start, peak, end, comp[peak[0]], 0, 0));
								
								//for(int k = comp_range[peak[0]][0]; k <= comp_range[peak[1]][1]; k++) valid[k] = true;
								
								if(comp[peak[0]] > 0.001)
								for(int k = comp_range[start][0]; k <= comp_range[end][1]; k++) valid[k] = true;
								//for(int k = max_idx[start]; k <= max_idx[end]; k++) valid[k] = true;
							}
						}
						peaks_list.add(peaks);
						valid_list.add(valid);
						*/


    					
    					
    					/*
    					double rms2 = 0;
    					double[] mag = new double[n_fft / 2];
    					double[] real_mag = new double[mag.length];
    					double[] freq = new double[mag.length];
    					double[] phase = new double[mag.length];
    					
    					
    					
    					//val[0] = complex[0].abs();
    					//val[n_fft / 2] = complex[n_fft / 2].abs();
    					//rms2 += val[0] * val[0];
    					for(int j = 0; j < n_fft / 2; j++) {
    						mag[j] = complex[j].abs() * eq[j];
    						real_mag[j] = complex[j].abs();
    						phase[j] = Math.atan2(complex[j].getImaginary(), complex[j].getReal());
    						freq[j] = j * bin_Hz;
    						if(phase[j] < 0) phase[j] += 2 * Math.PI;
    						rms2 += mag[j] * mag[j];
    					}
    					
    					
    					rms2 = Math.sqrt(rms2 / (n_fft / 2));
    					System.out.println(rms + ", " + rms2);
    					
    					for(int j = 0; j < mag.length; j++) if(mag_max < mag[j]) mag_max = mag[j];
    					mag_list.add(mag);
    					
    					//for(int j = 0; j < mag.length; j++) if(mag_max < real_mag[j]) mag_max = real_mag[j];
    					
    					//for(int j = 1; j < mag.length - 1; j++) if(mag[j] < Math.max(mag[j - 1], mag[j + 1])) real_mag[j] = 0;
    					
    					
    					//mag_list.add(real_mag);
    					
    					
    					freq_list.add(freq);
    					phase_list.add(phase);
    					*/


    					/*
    					ArrayList<Peak> peaks = new ArrayList<>();
						for(int j = 0; j < mag.length; j++) {
							if(mag[j] > 0) {
								int start = j;
								int[] peak = new int[2];
								peak[0] = j;
								while(j + 1 < mag.length && mag[j] <= mag[j + 1]) {
									if(mag[j] < mag[j + 1]) peak[0] = j + 1;
									j++;
								}
								peak[1] = j;
								while(j + 1 < mag.length && 0 < mag[j + 1] && mag[j + 1] <= mag[j]) {
									j++;
								}
								int end = j;
								
								//double value = dbfs[peak[0]];
								double v = 0;
								for(int k = start; k <= end; k++) v += real_mag[k];
								int mid = (peak[0] + peak[1]) / 2;
								
								double ph = Math.atan2(complex[mid].getImaginary(), complex[mid].getReal());
								double fr = mid * bin_Hz;
								//double value = Math.max(0, 20 * Math.log10(scal[peak[0]]) + dbNorm) / dbNorm;
								//if(valid[peak[0]])
								//if(end - start + 1 > mel_px)
								peaks.add(new Peak(frame, start, peak, end, v, ph, fr));
								//sum += value;
								//if(comp_max < value) comp_max = value;
							}
						}
						peaks_list.add(peaks);
						*/

    					/*
    					double[] mel = new double[mel_range.length];
    					double[] eq_mel = new double[mel_range.length];
    					double[] freq = new double[mel_range.length];
    					double[] phase = new double[mel_range.length];
    					for(int j = 0; j < mel_range.length; j++) {
    						for(int k = mel_range[j][0]; k <= mel_range[j][1]; k++) {
    							double eqmel = complex[k].abs() * eq[k];
    							if(eq_mel[j] < eqmel) {
    								freq[j] = k * bin_Hz;
    								mel[j] = complex[k].abs();
    								eq_mel[j] = eqmel;
    								phase[j] = Math.atan2(complex[k].getImaginary(), complex[k].getReal());
    							}
    						}
    						if(phase[j] < 0) phase[j] += 2 * Math.PI;
    						
    						
    						//if(mag_max < mel[j]) mag_max = mel[j];
    						
    						if(mag_max < eq_mel[j]) mag_max = eq_mel[j];
    					}
    					mag_list.add(eq_mel);
    					//mag_list.add(mel);
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
    							//freq[j] = k * bin_Hz;
    							//phase[j] = Math.atan2(complex[k].getImaginary(), complex[k].getReal());
    						}
    						if(mag_max < mel[j]) mag_max = mel[j];
    					}
    					mag_list.add(mel);
    					phase_list.add(phase);
    					freq_list.add(freq);
    					*/

    					
    					/*
    					double[] mag = new double[n_fft / 2];
    					double[] freq = new double[mag.length];
    					double[] phase = new double[mag.length];
    					for(int j = 1; j < mel.length - 1; j++) {
    						if(mel[j] > Math.max(mel[j - 1], mel[j + 1])) {
    							double fr = mel_freq[j];
    							int freq_idx = (int)(fr / bin_Hz);
    							mag[freq_idx] = complex[freq_idx].abs();
    							freq[freq_idx] = fr;
    							phase[freq_idx] = Math.atan2(complex[freq_idx].getImaginary(), complex[freq_idx].getReal());
    							
    							if(mag_max < mag[freq_idx]) mag_max = mag[freq_idx];
    						}
    					}
    					mag_list.add(mag);
    					phase_list.add(phase);
    					freq_list.add(freq);
    					*/
    					

    					/*
    					double[] mag = new double[n_fft / 2];
    					double max = 0;
    					for(int j = 0; j < mag.length; j++) {
    						mag[j] = complex[j].abs();
    						if(max < mag[j]) max = mag[j];
    					}
    					*/
    					/*
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
								if(value[peak[0]] > rms2) peaks.add(new Peak(frame, start, peak, end, value[peak[0]]));
								//sum += value;
								//if(comp_max < value) comp_max = value;
							}
						}
						peaks_list.add(peaks);
						*/
    					
    					
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
                        
                        if(frame == 1500) bb = true;
                    }
            	}
            	if(bb) break;
            }
            System.out.println(System.currentTimeMillis() - stamp);
            
            
            BufferedImage peak = new BufferedImage(peaks_list.size(), comp_range.length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < peaks_list.size(); i++) {
            	ArrayList<Peak> peaks = peaks_list.get(i);
            	for(int j = 0; j < peaks.size(); j++) {
            		for(int k = peaks.get(j).start; k <= peaks.get(j).end; k++) {
            			int col = 0;
            			int degree = (int)(peaks.get(j).value / mag_eq_max * 0xFF);
            			//degree = 0xFF;
            			if(peaks.get(j).peak[0] <= k && k <= peaks.get(j).peak[1]) col = (degree << 24 | degree << 16);
            			else if(k < peaks.get(j).peak[0]) col = (degree << 24 | degree << 8);
            			else if(k > peaks.get(j).peak[1]) col = (degree << 24 | degree);

            			peak.setRGB(i, comp_range.length - 1 - k, col);
            		}
            	}
            }
            ImageIO.write(peak, "PNG", new File("output_peak.png"));//
            
            /*
            BufferedImage peak = new BufferedImage(peaks_list.size(), mag_pq_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < peaks_list.size(); i++) {
            	ArrayList<Peak> peaks = peaks_list.get(i);
            	for(int j = 0; j < peaks.size(); j++) {
            		for(int k = peaks.get(j).start; k <= peaks.get(j).end; k++) {
            			int col = 0;
            			int degree = (int)(peaks.get(j).value / mag_eq_max * 0xFF);
            			//degree = 0xFF;
            			if(peaks.get(j).peak[0] <= k && k <= peaks.get(j).peak[1]) col = (degree << 24 | degree << 16);
            			else if(k < peaks.get(j).peak[0]) col = (degree << 24 | degree << 8);
            			else if(k > peaks.get(j).peak[1]) col = (degree << 24 | degree);

            			peak.setRGB(i, mag_pq_list.get(0).length - 1 - k, col);
            		}
            	}
            }
            ImageIO.write(peak, "PNG", new File("output_peak.png"));//
            */
            
            BufferedImage energy = new BufferedImage(peak_sum_list.size(), 500, BufferedImage.TYPE_INT_RGB);
            for(int i = 1; i < peak_sum_list.size() - 1; i++) {
            	int height = (int)((energy.getHeight() - 1) * (Math.abs(peak_sum_list.get(i)) / max_peak_sum));
            	
            	//boolean peak = smooth_energy_list.get(i) > Math.max(smooth_energy_list.get(i - 1), smooth_energy_list.get(i + 1));
            	
            	for(int j = 0; j <= height; j++) {
            		//energy.setRGB(i, (energy.getHeight() - 1) - j, peak ? 0xFFFF0000 : 0xFFFFFFFF);
            		energy.setRGB(i, (energy.getHeight() - 1) - j, peak_sum_list.get(i) < 0 ? 0xFFFF0000 : 0xFFFFFFFF);
            	}
            }
            ImageIO.write(energy, "PNG", new File("output_energy_list.png"));
            
            /*
            BufferedImage valid_img = new BufferedImage(valid_list.size(), valid_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < valid_list.size(); i++) {
            	boolean[] valid = valid_list.get(i);
            	
            	for(int j = 1; j < valid.length - 1; j++) {
            		//int col = (int)(valid[j] / mag_max * 0xFF);
            		//dbfs.setRGB(i, mag_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            		if(valid[j]) valid_img.setRGB(i, valid_list.get(0).length - 1 - j, 0xFFFFFFFF);
            	}
            }
            ImageIO.write(valid_img, "PNG", new File("output_valid.png"));//
            */

            /*
            BufferedImage wave_img = new BufferedImage(wave_list.size(), wave_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < wave_list.size(); i++) {
            	Wave[] wave = wave_list.get(i);
            	
            	for(int j = 0; j < wave.length; j++) {
            		int col = (int)(wave[j].amplitude / mag_max * 0xFF);
            		wave_img.setRGB(i, wave_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            	}
            }
            ImageIO.write(wave_img, "PNG", new File("output_wave.png"));
            */        
            
            /*
            BufferedImage dbfs = new BufferedImage(mag_list.size(), mag_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mag_list.size(); i++) {
            	double[] mag = mag_list.get(i);
            	
            	for(int j = 1; j < mag.length - 1; j++) {
            		int col = (int)(mag[j] / mag_max * 0xFF);
            		//if(col > 0) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		if(col > 0) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFFFFFF);
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            */

            /*
            BufferedImage dbfs4 = new BufferedImage(mel_list.size(), mel_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mel_list.size(); i++) {
            	double[] mel = mel_list.get(i);
            	
            	for(int j = 1; j < mel.length - 1; j++) {
            		int col = (int)(mel[j] / mel_max * 0xFF);
            		dbfs4.setRGB(i, mel_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs4, "PNG", new File("output_dbfs4.png"));//
            */
            
            BufferedImage phase = new BufferedImage(phase_list.size(), phase_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < phase_list.size(); i++) {
            	double[] phase_arr = phase_list.get(i);
            	
            	for(int j = 0; j < phase_arr.length; j++) {
            		int col = (int)(phase_arr[j] / (2 * Math.PI) * 0xFF);
            		phase.setRGB(i, phase_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(phase, "PNG", new File("output_phase.png"));//

            BufferedImage dbfs2 = new BufferedImage(mag_eq_list.size(), mag_eq_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mag_eq_list.size(); i++) {
            	double[] mag_eq = mag_eq_list.get(i);
            	
            	for(int j = 1; j < mag_eq.length - 1; j++) {
            		int col = (int)(mag_eq[j] / mag_eq_max * 0xFF);
            		dbfs2.setRGB(i, mag_eq_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs2, "PNG", new File("output_dbfs2.png"));//

            /*
            BufferedImage dbfs3 = new BufferedImage(mag_pq_list.size(), mag_pq_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mag_pq_list.size(); i++) {
            	double[] mag_pq = mag_pq_list.get(i);
            	
            	for(int j = 1; j < mag_pq.length - 1; j++) {
            		int col = (int)(mag_pq[j] / mag_pq_max * 0xFF);
            		dbfs3.setRGB(i, mag_pq_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs3, "PNG", new File("output_dbfs3.png"));//
            */
            
            /*
            double[][] x_ = new double[][] {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
            double[][] y_ = new double[][] {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
            
            double[][] sob = new double[mag_eq_list.size()][mag_eq_list.get(0).length];
            double sob_max = 0;
            
            for(int i = 1; i < mag_eq_list.size() - 1; i++) {
            	for(int j = 1; j < mag_eq_list.get(0).length  - 1; j++) {
            		double sobel_x = 0;
            		double sobel_y = 0;
            		for(int x = -1; x <= 1; x++) {
            			for(int y = -1; y <= 1; y++) {
            				sobel_x += mag_eq_list.get(i + x)[j + y] * x_[x + 1][y + 1];
            				sobel_y += mag_eq_list.get(i + x)[j + y] * y_[x + 1][y + 1];
            			}
            			
            			double rad = Math.atan2(sobel_y, sobel_x);
		        		if(rad < 0) rad += Math.PI * 2;
		        		//sob[i][j] = rad;
		        		sob[i][j] = Math.sqrt(sobel_x * sobel_x + sobel_y * sobel_y);
		        		if(sob_max < sob[i][j]) sob_max = sob[i][j];
            		}
            	}
            }
            
            BufferedImage sob_img = new BufferedImage(sob.length, sob[0].length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < sob.length; i++) {
            	for(int j = 0; j < sob[0].length; j++) {
            		int col = (int)(sob[i][j] / sob_max * 0xFF);
            		//if(col > 0) sob_img.setRGB(i, sob[0].length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		if(col > 1) sob_img.setRGB(i, sob[0].length - 1 - j, 0xFFFFFFFF);
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(sob_img, "PNG", new File("output_sob.png"));//
            */
            

            BufferedImage gradient_img = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < gradient_list.size(); i++) {
            	double[][] grad = gradient_list.get(i);
            	for(int j = 0; j < grad.length; j++) {
            		int col = (int)(grad[j][1] / gradient_max * 0xFF);
            		gradient_img.setRGB(i, grad.length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(gradient_img, "PNG", new File("output_gradient.png"));//
            
            

            BufferedImage soh = new BufferedImage(sum_of_harmonic_list.size(), sum_of_harmonic_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < soh.getWidth(); i++) {
            	for(int j = 0; j < soh.getHeight(); j++) {
        			int col = (int)(sum_of_harmonic_list.get(i)[j] / max_sum_of_harmonic * 0xFF);
        			soh.setRGB(i, sum_of_harmonic_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
        			if(max_idx_of_harmonic_list.get(i) == j) soh.setRGB(i, sum_of_harmonic_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(soh, "PNG", new File("output_soh.png"));//

            /*
            BufferedImage ang = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            BufferedImage deg = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            
            double a = 90. / 180. * Math.PI; //160.
            for(int i = 0; i < ang.getWidth(); i++) {
            	for(int j = 1; j < gradient_list.get(0).length - 1; j++) {
            		//int degree = 0xFF;
            		int degree = Math.min(0xFF, (int)(gradient_list.get(i)[j][1] * 0xFF));
            		//if(gradient_list.get(i)[j][1] > 1. / 65536) degree = 0xFF;
            		
            		//degree = 0xFF;
            		
        			int col = 0;
        			
        			//int sec = (int)Math.floor((gradient_list.get(i)[j][0] + Math.PI / 4) / (Math.PI / 2)) % 4;
        			
        			double angle = gradient_list.get(i)[j][0] + a / 2;
        			int sec = (int)(Math.floor(angle / Math.PI) * 2 + (Math.floor(angle % Math.PI / a) > 0 ? 1 : 0)) % 4;
        			
        			if(sec == 0) col = (degree << 24 | degree << 16);
            		else if(sec == 1) col = (degree << 24 | degree << 8);
            		else if(sec == 2) col = (degree << 24 | degree);
            		else if(sec == 3) col = (degree << 24 | degree << 16 | degree);
        			
        			
        			//if(gradient_list.get(i)[j][1] > 10) ang.setRGB(i, gradient_list.get(0).length - 1 - j, col);
        			if(gradient_list.get(i)[j][1] > 0.) ang.setRGB(i, gradient_list.get(0).length - 1 - j, col);
        			
        			//if(i < valid_list.size() && valid_list.get(i)[j]) ang.setRGB(i, gradient_list.get(0).length - 1 - j, col);
        			
        			int d = Math.min(0xFF, (int)(gradient_list.get(i)[j][1] * 0xFF));
        			
        			//if(gradient_list.get(i)[j][1] > 0) System.out.println(gradient_list.get(i)[j][1]);
        			//if(d > 10) d = 0xFF;
        			deg.setRGB(i, gradient_list.get(0).length - 1 - j, (d << 24 | d << 16 | d << 8 | d));
            	}
            }
            ImageIO.write(ang, "PNG", new File("output_gradient_angle.png"));//
            ImageIO.write(deg, "PNG", new File("output_gradient_magnitude.png"));//
            */


            
            
            /*
            int frame_samples = (int)(frame_sec * sampleRate);
            double sampleMax = 0;
            double[] samplesArray = new double[frame_samples * peaks_list.size()];
            for(int i = 0; i < peaks_list.size(); i++) {
            	ArrayList<Peak> peaks = peaks_list.get(i);
            	for(int j = 0; j < frame_samples; j++) {
            		double common = 2 * Math.PI * (j + 0.) / sampleRate;
            		for(int k = 0; k < peaks.size(); k++) samplesArray[i * frame_samples + j] += Math.sin(common * peaks.get(k).freq + peaks.get(k).phase) * peaks.get(k).value;
            		sampleMax = Math.max(sampleMax, Math.abs(samplesArray[i * frame_samples + j]));
            	}
            	System.out.println(i);
            }
            byte[] pcmBytesArray = new byte[samplesArray.length * 2];
            for(int i = 0; i < samplesArray.length; i++) {
            	pcmBytesArray[i * 2] = (byte)((int)(samplesArray[i] / sampleMax * 32767) & 0xFF);
            	pcmBytesArray[i * 2 + 1] = (byte)((int)(samplesArray[i] / sampleMax * 32767) >> 8 & 0xFF);
            }
            */
            
            /*
            int frame_samples = (int)(frame_sec * sampleRate);
            double sampleMax = 0;
            double[] samplesArray = new double[frame_samples * mag_list.size()];
            for(int i = 0; i < mag_list.size() - 1; i++) {
            	double[] mag = mag_list.get(i);
            	double[] next_mag = mag_list.get(i + 1);
            	for(int j = 0; j < frame_samples; j++) {
            		double common = 2 * Math.PI * (j + 0.) / sampleRate;
            		for(int k = 0; k < mag.length; k++) {
            			double mag_lerp = mag[k] + (next_mag[k] - mag[k]) * ((double)j / frame_samples);
            			samplesArray[i * frame_samples + j] += Math.cos(common * comp_freq[k]) * mag_lerp;
            		}
            		sampleMax = Math.max(sampleMax, Math.abs(samplesArray[i * frame_samples + j]));
            	}
            	System.out.println(i);
            }
            */
            
 
            int frame_samples = (int)(frame_sec * sampleRate);
            double sampleMax = 0;
            double[] samplesArray = new double[frame_samples * mag_list.size()];
            double[] angle = new double[overtone_count];
            for(int i = 0; i < mag_list.size() - 1; i++) {
            	double[] mag = mag_list.get(i);
            	int max_idx = max_idx_of_harmonic_list.get(i);
            	
            	double[] mag2 = mag_list.get(i + 1);
            	int max_idx2 = max_idx_of_harmonic_list.get(i + 1);
            	//int[] max_idx2 = max_list.get(i);
            	
            	
            	for(int j = 0; j < frame_samples; j++) {
            		for(int k = 0; k < angle.length; k++) {
            			int l = max_idx + (int)Math.round(12 * mel_comp * Math.log(1 + k * overtone_coef) / Math.log(2));
            			
            			/*
            			if(l < comp_freq.length) {
            				angle[k] += 2 * Math.PI * comp_freq[l] / sampleRate;
	                		samplesArray[i * frame_samples + j] += Math.cos(angle[k]) * mag[l];
            			}
            			*/
            			
            			int l2 = max_idx2 + (int)Math.round(12 * mel_comp * Math.log(1 + k * overtone_coef) / Math.log(2));
						if(l < comp_freq.length && l2 < comp_freq.length) {
							//angle[k] += 2 * Math.PI * (max_idx2[l] * bin_Hz) / sampleRate;
							angle[k] += 2 * Math.PI * (comp_freq[l] + (comp_freq[l2] - comp_freq[l]) * ((double)j / frame_samples)) / sampleRate;
	                		samplesArray[i * frame_samples + j] += Math.cos(angle[k]) * (mag[l] + (mag2[l2] - mag[l]) * ((double)j / frame_samples));
						}
            		}
            		sampleMax = Math.max(sampleMax, Math.abs(samplesArray[i * frame_samples + j]));
            	}
            	//max_idx는 0배음으로 칠것임
            }

            byte[] pcmBytesArray = new byte[samplesArray.length * 2];
            for(int i = 0; i < samplesArray.length; i++) {
            	pcmBytesArray[i * 2] = (byte)((int)(samplesArray[i] / sampleMax * 32767) & 0xFF);
            	pcmBytesArray[i * 2 + 1] = (byte)((int)(samplesArray[i] / sampleMax * 32767) >> 8 & 0xFF);
            }


            
            pcmToWav((int)sampleRate, channelCount, 16, pcmBytesArray, "output.wav");

            
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
}
