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
            double win_sec = 0.08;
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


            
            /*
            double[] hz = new double[] { 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500 };
            double[] af = new double[] { 0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349, 0.330, 0.315, 0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250, 0.246, 0.244, 0.243, 0.243, 0.243, 0.242, 0.242, 0.245, 0.254, 0.271, 0.301 };
            double[] lu = new double[] { -31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3, -8.1, -6.2, -4.5, -3.1, -2.0, -1.1, -0.4, 0.0, 0.3, 0.5, 0.0, -2.7, -4.1, -1.0, 1.7, 2.5, 1.2, -2.1, -7.1, -11.2, -10.7, -3.1 };
            double[] tf = new double[] { 78.5, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9, 14.4, 11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3, -4.2, -6.0, -5.4, -1.5, 6.0, 12.6, 13.9, 12.3 };
            double[] eq = new double[n_fft / 2];
            for(int i = 0; i < eq.length; i++) {
            	double freq = i * sampleRate / n_fft;
            	int low = 0;
            	int high = hz.length;
            	while(low < high) {
            		int mid = (low + high) / 2;
            		if(freq < hz[mid]) high = mid;
            		else low = mid + 1;
            	}
            	if(0 < low && low < hz.length) {
            		double lerp = (freq - hz[low - 1]) / (hz[low] - hz[low - 1]);
            		double new_af = af[low - 1] + (af[low] - af[low - 1]) * lerp;
            		double new_tf = tf[low - 1] + (tf[low] - tf[low - 1]) * lerp;
            		double new_lu = lu[low - 1] + (lu[low] - lu[low - 1]) * lerp;
            		double Ln = 40;
            		double Af = 0.00447 * (Math.pow(10, 0.025 * Ln) - 1.15) + Math.pow(0.4 * Math.pow(10, (new_tf + new_lu) / 10 - 9), new_af);
            		double Lp = 10 / new_af * Math.log10(Af) - new_lu + 94;
            		
            		//eq[i] = Lp / 20;
            		eq[i] = Math.pow(10, -Lp / 20);
            		
            		System.out.println((int)(freq) + " : " + eq[i]);
            		if(equal_loudness_max < eq[i]) equal_loudness_max = eq[i];
            	}
            }
            */

            System.out.println(equal_loudness_max);

            for(int i = 0; i < eq.length; i++) eq[i] /= equal_loudness_max;
            BufferedImage img = new BufferedImage(eq.length, 500, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < img.getWidth(); i++) {
            	img.setRGB(i, (int)((img.getHeight() - 1) * (1 - eq[i])), 0xFFFFFFFF);
            }
            ImageIO.write(img, "PNG", new File("output_eq.png"));//
            

            /*
            int mel_comp = 1;
            int[][] mel_range2 = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp)][2];
            //int[][] mel_range2 = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp)][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key = (int)Math.round(Math.max(-1, 12 * Math.log((double)i * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	if(0 <= key && key < mel_range2.length) {
            		if(mel_range2[key][0] == 0) mel_range2[key][0] = i;
            		mel_range2[key][1] = i;
            	}
            }
            */
            
            int mel_comp = 10;
            int[][] mel_range = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp) + 1][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key1 = (int)Math.round(Math.max(-1, 12 * Math.log((i + 0.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	int key2 = (int)Math.round(Math.max(-1, 12 * Math.log((i + 1.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	for(int j = Math.max(0, key1); j <= key2; j++) {
            		if(mel_range[j][0] == 0) mel_range[j][0] = i;
            		mel_range[j][1] = i;
            	}
            }

            
            
            double bin_Hz = sampleRate / n_fft;

            
            /*
            int mel_comp = 10;
            int keyCount = 100;
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

            
            ArrayList<ArrayList<Peak>> peaks_list = new ArrayList<>();
            
            FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            
            List<Double> frame_list = new ArrayList<>();
            
            List<List<Mel>> mel_frames = new ArrayList<>();

            
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
                    	
                    	
    					for(int j = 0; j < n_fft; j++) {
    						complex[j] = new Complex(j < sampleBuffer.length ? sampleBuffer[j] * Math.sqrt(win[j]) : 0, 0);
    					}
    					complex = fft.transform(complex, TransformType.FORWARD);
    					
    					
    					double rms = 0;
    					for(int j = 0; j < win.length; j++) {
    						rms += sampleBuffer[j] * sampleBuffer[j] * win[j];
    					}
    					rms = Math.sqrt(rms);

    					/*
    					double rms2 = 0;
    					for(int j = 0; j < n_fft; j++) {
    						rms2 += complex[j].abs() * complex[j].abs();
    					}
    					rms2 /= n_fft;
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
                        
                        if(frame == 1000) bb = true;
                    }
            	}
            	if(bb) break;
            }
            System.out.println(System.currentTimeMillis() - stamp);
            
            /*
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
            */
            
            BufferedImage dbfs = new BufferedImage(mag_list.size(), mag_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mag_list.size(); i++) {
            	double[] mag = mag_list.get(i);
            	
            	/*
            	//if(i == 63) {
            	if(i == 325) {
            		int[] max_mel = max_list.get(i);
            		double[] phase = phase_list.get(i);
            		double[] angle = new double[phase.length];
            		double[] freq = freq_list.get(i);
            		
            		BufferedImage im = new BufferedImage(mag.length, 500, BufferedImage.TYPE_INT_RGB);
            		double max = 0;
            		for(int j = 0; j < mag.length; j++) if(max < mag[j]) max = mag[j];
            		
            		for(int j = 0; j < phase.length; j++) {
            			angle[j] = phase[j] / Math.PI * 180;
            			//if(j == 145 || j == 294 || j == 442) System.out.println(j + " : " + freq[j] / freq[]);
            			//if(j == 144 || j == 293) System.out.println(j + " : " + angle[j]);
            			//if(j == 146 || j == 295) System.out.println(j + " : " + angle[j]);
            		}
            		
            		System.out.println(freq[145] / freq[145]);
            		System.out.println(freq[294] / freq[145]);
            		System.out.println(freq[442] / freq[145]);
            		
            		for(int j = 0; j < mag.length; j++) {
            			int height = (int)(mag[j] / max * (im.getHeight() - 1));
            			for(int k = 0; k <= height; k++) im.setRGB(j, (im.getHeight() - 1) - k, 0xFFFFFFFF);
            		}
            		
            		for(int j = 0; j < max_mel.length; j++) {
            			int max_idx = max_mel[j];
            			int height = (int)(mag[max_idx] / max * (im.getHeight() - 1));
            			for(int k = 0; k <= height; k++) im.setRGB(max_idx, (im.getHeight() - 1) - k, 0xFFFF0000);
            		}
            		
            		ImageIO.write(im, "PNG", new File("output_im" + n_fft  + ".png"));//
            	}
            	*/
            	
            	
            	for(int j = 1; j < mag.length - 1; j++) {
            		int col = (int)(mag[j] / mag_max * 0xFF);
            		dbfs.setRGB(i, mag_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            
            
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
            
            int frame_samples = (int)(frame_sec * sampleRate);
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
            		
            		for(int k = 0; k < mag.length; k++) samplesArray[i * frame_samples + j] += Math.cos(common * freq[k] + phase[k]) * mag[k];
            		
            		//for(int k = 0; k < mag.length; k++) samplesArray[i * frame_samples + j] += mag[k] * Math.cos(common * freq[k]) * Math.cos(phase[k]) - mag[k] * Math.sin(common * freq[k]) * Math.sin(phase[k]);
            		
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
