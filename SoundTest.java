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
			FileInputStream fis = new FileInputStream(new File("sample/sample1.pcm"));
			int sampleRateOrigin = 44100;
			double sampleRate = sampleRateOrigin;

            double frame_sec = 0.02;
            //System.out.println(frame_sec);
            double win_sec = frame_sec * 4;
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
            //n_fft = 32768;
            
            System.out.println("1px : " + (sampleRate / n_fft) + "Hz");
            
            System.out.println("win : " + rect_win.length + ", n_fft : " + n_fft);
            
            
            double[] sampleBuffer = new double[rect_win.length];
            int sampleOffset = sampleBuffer.length - frame_len;
            

            int mel_comp = 10;
            int[][] comp_range = new int[(int)Math.round(12 * Math.log(sampleRate * 0.5 / A0) / Math.log(2) * mel_comp) + 1][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key1 = (int)Math.round(Math.max(-1, 12 * Math.log((i + 0.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	int key2 = (int)Math.round(Math.max(-1, 12 * Math.log((i + 1.) * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp);
            	for(int j = Math.max(0, key1); j <= key2; j++) {
            		if(comp_range[j][0] == 0) comp_range[j][0] = i;
            		comp_range[j][1] = i;
            	}
            }
            
            double bin_Hz = sampleRate / n_fft;
            
            double equal_loudness_max = 0;
            double[] eq = new double[n_fft / 2];
            for(int i = 0; i < eq.length; i++) {
				eq[i] = equal_loudness(i * bin_Hz);
            	if(equal_loudness_max < eq[i]) equal_loudness_max = eq[i];
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

            ArrayList<double[]> mel_list = new ArrayList<>();
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
            
            List<Wave[]> wave_list = new ArrayList<>();
            
            
            List<boolean[]> valid_list = new ArrayList<>();
            
            List<List<Wave>> mel_frames = new ArrayList<>();

            
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
    					for(int j = 0; j < n_fft; j++) complex[j] = new Complex(j < sampleBuffer.length ? sampleBuffer[j] * Math.sqrt(rect_win[j]) : 0, 0);
    					complex = fft.transform(complex, TransformType.FORWARD);
    					
    					//Wave[] wave = new Wave[n_fft / 2];
    					//for(int j = 0; j < wave.length; j++) wave[j] = new Wave(j * bin_Hz, complex[j].abs(), Math.atan2(complex[j].getImaginary(), complex[j].getReal()));
    					//wave_list.add(wave);
    					
    					double[] mag = new double[n_fft / 2];
    					for(int j = 0; j < mag.length; j++) mag[j] = complex[j].abs();
    					
    					
    					sobel.Add(frame, mag, false);
    					while(sobel.rslt_list.size() > 0) {
    						double[][] gradient = sobel.rslt_list.remove(0);
    						gradient_list.add(gradient);
    					}
    					
    					
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
    					
    					
    					double[] sum_of_harmonic = new double[mel.length];
    					
    					int max_idx_of_harmonic = 0;
    					
						for(int k = 0; k < mel.length; k++) {
							for(int j = 0; j < overtone_count; j++) {
		                		int l = k + (int)Math.round(12 * mel_comp * Math.log(1 + j * 1.0) / Math.log(2));
		                		
								if(l < mel.length) sum_of_harmonic[k] += mel[l]; //comp[l];
		                    }
							if(sum_of_harmonic[max_idx_of_harmonic] < sum_of_harmonic[k]) max_idx_of_harmonic = k;
						}
						if(max_sum_of_harmonic < sum_of_harmonic[max_idx_of_harmonic]) max_sum_of_harmonic = sum_of_harmonic[max_idx_of_harmonic];
						
						
						max_idx_of_harmonic_list.add(max_idx_of_harmonic);
						sum_of_harmonic_list.add(sum_of_harmonic);
    					
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
                        
                        //if(frame == 10) bb = true;
                    }
            	}
            	if(bb) break;
            }
            System.out.println(System.currentTimeMillis() - stamp);
            

            /*
            BufferedImage peak = new BufferedImage(peaks_list.size(), comp_range.length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < peaks_list.size(); i++) {
            	ArrayList<Peak> peaks = peaks_list.get(i);
            	Wave[] wave = wave_list.get(i);
            	for(int j = 0; j < peaks.size(); j++) {
            		for(int k = peaks.get(j).start; k <= peaks.get(j).end; k++) {
            			int col = 0;
            			int degree = (int)(peaks.get(j).value / mag_max * 0xFF);
            			
            			if(peaks.get(j).peak[0] <= k && k <= peaks.get(j).peak[1]) col = (degree << 24 | degree << 16);
            			else if(k < peaks.get(j).peak[0]) col = (degree << 24 | degree << 8);
            			else if(k > peaks.get(j).peak[1]) col = (degree << 24 | degree);

            			peak.setRGB(i, comp_range.length - 1 - k, col);
            		}
            	}
            }
            ImageIO.write(peak, "PNG", new File("output_peak.png"));//
            */
            
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
            
            
            
            BufferedImage dbfs = new BufferedImage(mag_list.size(), mag_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mag_list.size(); i++) {
            	double[] mag = mag_list.get(i);
            	
            	for(int j = 1; j < mag.length - 1; j++) {
            		int col = (int)(mag[j] / mag_max * 0xFF);
            		dbfs.setRGB(i, mag_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) dbfs.setRGB(i, mag_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            
            
            BufferedImage soh = new BufferedImage(sum_of_harmonic_list.size(), sum_of_harmonic_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < soh.getWidth(); i++) {
            	for(int j = 0; j < soh.getHeight(); j++) {
        			int col = (int)(sum_of_harmonic_list.get(i)[j] / max_sum_of_harmonic * 0xFF);
        			soh.setRGB(i, sum_of_harmonic_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
        			if(max_idx_of_harmonic_list.get(i) == j) soh.setRGB(i, sum_of_harmonic_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(soh, "PNG", new File("output_soh.png"));//
            
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
            double sampleMax = 0;
            double[] samplesArray = new double[frame_samples * mag_list.size()];
            double[] angle = new double[overtone_count];
            for(int i = 0; i < mag_list.size(); i++) {
            	double[] mag = mag_list.get(i);
            	int max_idx = max_idx_of_harmonic_list.get(i);
            	
            	int[] max_idx2 = max_list.get(i);
            	
            	
            	for(int j = 0; j < frame_samples; j++) {
            		for(int k = 0; k < angle.length; k++) {
            			int l = max_idx + (int)Math.round(12 * mel_comp * Math.log(1 + k * 1.0) / Math.log(2));
						if(l < max_idx2.length) {
							
							angle[k] += 2 * Math.PI * (max_idx2[l] * bin_Hz) / sampleRate;
							
	            			//angle[k] += 2 * Math.PI * mel_freq[l] / sampleRate;
	                		samplesArray[i * frame_samples + j] += Math.cos(angle[k]) * (mag[max_idx2[l]]);
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
