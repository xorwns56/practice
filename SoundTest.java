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
			
			System.out.println(i);
		}
		
		return pcmBytesArray;
	}
	
	static final double A0 = 13.75;
	
	public static double i0(double x){
		return (Math.cosh(x) + 2 * (Math.cosh(0.970941817426052 * x) + Math.cosh(0.8854560256532099 * x) + Math.cosh(0.7485107481711011 * x) + Math.cosh(0.5680647467311558 * x) + Math.cosh(0.3546048870425356 * x) + Math.cosh(0.120536680255323 * x))) / 13;
	}
	
	public static double[] kaiserWindowedSinc(double fp, double fs, double A) {
		double b = fs - fp;
		double fc = fp + b * 0.5;
		int N = (int)((A - 8) / (2.285 * Math.PI * b)) + 1;
		if(N % 2 == 0) N++;
		double beta = 0;
		if(A > 50) beta = 0.1102 * (A - 8.7);
		else if(21 <= A && A <= 50) beta = 0.5842 * Math.pow(A - 21, 0.4) + 0.07886 * (A - 21);
		double[] window = new double[N];
		double window_sum = 0;
		for(int i = 0; i < N; i++) {
			window[i] = i0(beta * Math.sqrt(1 - Math.pow((2. * i / (N - 1) - 1), 2))) / i0(beta);
        	double x = 2 * Math.PI * fc * (i - N / 2);
        	window[i] *= x != 0 ? Math.sin(x) / x : 1;
        	window_sum += window[i];
		}
		for(int i = 0; i < N; i++) window[i] /= window_sum;
		return window;
    }
	
	public static double[] blur2(double[] src, int r) {
		double[] tmp = new double[src.length];
		double div = 1. / (2 * r + 1);
		double value = (r + 1) * src[0];
		for(int j = 0; j < r; j++) value += src[j];
		for(int j = 0; j < r + 1; j++) {
			value += src[j + r] - src[0];
			tmp[j] = value * div;
		}
		for(int j = r + 1; j < src.length - r; j++) {
			value += src[j + r] - src[j - (r + 1)];
			tmp[j] = value * div;
		}
		for(int j = src.length - r; j < src.length; j++) {
			value += src[src.length - 1] - src[j - (r + 1)];
			tmp[j] = value * div;
		}
		return tmp;
	}
	
	private static double[] make_curve(double sigma) {
		double sigma2 = 2 * sigma * sigma;
		int n = (int)Math.ceil(Math.sqrt(-sigma2 * Math.log(1. / 255))) * 2;
		if(n % 2 == 0) n++;
		int r = n / 2;
		double[] curve = new double[n];
		curve[r] = 255;
		double total = curve[r];
		for(int i = 1; i <= r; i++) {
			curve[r - i] = curve[r + i] = (int)(Math.exp(-(i * i) / sigma2) * 255);
			total += curve[r - i] * 2;
		}
		for(int i = 0; i < curve.length; i++) curve[i] /= total;
		return curve;
	}
	
	private static double[] gaus(double[] value, double[] kernel) {
		double[] tmp = new double[value.length];
		for(int i = 0; i < value.length; i++) {
			double sum = 0;
			int start = i - kernel.length / 2;
			int end = i + kernel.length / 2;
			for(int j = start; j < 0; j++) sum += value[0] * kernel[j - start];
			if(end >= value.length) {
				for(int j = value.length; j <= end; j++) sum += value[value.length - 1] * kernel[j - start];
				end = value.length - 1;
			}
			for(int j = Math.max(0, start); j <= end; j++) sum += value[j] * kernel[j - start];
			tmp[i] = sum;
		}
		return tmp;
	}
	
	private static class IndexValuePair{
		public int index;
		public double value;
		public IndexValuePair(int index, double value) {
			this.index = index;
			this.value = value;
		}
	}
	
	public static double[] median(double[] value, int r) {
		double[] tmp = new double[value.length];
		ArrayList<IndexValuePair> list = new ArrayList<>();
		for(int j = -r; j < 0; j++) list.add(new IndexValuePair(j, value[0]));
		for(int j = 0; j < value.length + r; j++) {
			IndexValuePair data = new IndexValuePair(j, value[Math.min(value.length - 1, j)]);
			int low = 0;
			int high = list.size();
			while(low < high) {
				int mid = (low + high) / 2;
				if((double)data.value < (double)list.get(mid).value) high = mid;
				else low = mid + 1;
			}
			list.add(low, data);
			if(j < r) continue;
			tmp[j - r] = (double)list.get(list.size() / 2).value;
			Iterator<IndexValuePair> iter = list.iterator();
			while(iter.hasNext()) {
				if(iter.next().index == j - r * 2) {
					iter.remove();
					break;
				}
			}
		}
		return tmp;
	}
	
	public static double[] blur(double[] src, int r) {
		double max_weight = 1. / (r + 1);
		double min_weight = max_weight / (r + 1);
		double[] tmp = new double[src.length];
		double moving_sum = 0;
		double weighted_moving_sum = 0;
		LinkedList<Double> list = new LinkedList<>();
		for(int j = 1; j <= r; j++) {
			weighted_moving_sum += src[0] * min_weight * j;
			list.add(src[0]);
			moving_sum += src[0];
		}
		for(int j = 0; j < src.length + r; j++) {
			double value = src[Math.min(src.length - 1, j)];
			weighted_moving_sum += value * max_weight;
			if(j < src.length) tmp[j] = weighted_moving_sum;
			list.add(value);
			moving_sum += value;
			weighted_moving_sum -= moving_sum * min_weight;
			moving_sum -= list.removeFirst();
			if(0 <= j - r) tmp[j - r] += moving_sum * max_weight - weighted_moving_sum;
		}
		return tmp;
	}

	
	public static class GausFilter{
		int v_r, h_r;
		double h_max_weight, h_min_weight;
        double[] h_moving_sum;
        double[] h_weighted_moving_sum;
        LinkedList<double[]> h_moving_window;
        LinkedList<double[]> h_rslt;
        LinkedList<double[]> rslt;
        public GausFilter(int v_r, int h_r) {
        	this.v_r = v_r;
        	this.h_r = h_r;
        	h_max_weight = 1. / (h_r + 1);
    		h_min_weight = h_max_weight / (h_r + 1);
    		h_moving_window = new LinkedList<>();
        	h_rslt = new LinkedList<>();
        	rslt = new LinkedList<>();
        }
        public void Add(int frame, double[] src, boolean last) {
        	double v_max_weight = 1. / (v_r + 1);
    		double v_min_weight = v_max_weight / (v_r + 1);
    		double[] v_rslt = new double[src.length];
    		double v_moving_sum = 0;
    		double v_weighted_moving_sum = 0;
    		LinkedList<Double> v_moving_window = new LinkedList<>();
    		for(int i = 1; i <= v_r; i++) {
    			v_weighted_moving_sum += src[0] * v_min_weight * i;
    			v_moving_window.add(src[0]);
    			v_moving_sum += src[0];
    		}
    		for(int i = 0; i < src.length + v_r; i++) {
    			double value = src[Math.min(src.length - 1, i)];
    			v_weighted_moving_sum += value * v_max_weight;
    			if(i < src.length) v_rslt[i] = v_weighted_moving_sum;
    			v_moving_window.add(value);
    			v_moving_sum += value;
    			v_weighted_moving_sum -= v_moving_sum * v_min_weight;
    			v_moving_sum -= v_moving_window.removeFirst();
    			if(0 <= i - v_r) v_rslt[i - v_r] += v_moving_sum * v_max_weight - v_weighted_moving_sum;
    		}
    		if(frame == 0) {
        		h_moving_sum = new double[src.length];
        		h_weighted_moving_sum = new double[src.length];
        		for(int i = 1; i <= h_r; i++) {
        			for(int j = 0; j < src.length; j++) {
        				h_weighted_moving_sum[j] += v_rslt[j] * h_min_weight * i;
        				h_moving_sum[j] += v_rslt[j];
        			}
        			h_moving_window.add(v_rslt);
        		}
        	}
        	for(int i = frame; i <= frame + (last ? h_r : 0); i++) {
        		for(int j = 0; j < src.length; j++) h_weighted_moving_sum[j] += v_rslt[j] * h_max_weight;
        		if(i == frame) h_rslt.add(h_weighted_moving_sum.clone());
        		h_moving_window.add(v_rslt);
        		for(int j = 0; j < src.length; j++) {
        			h_moving_sum[j] += v_rslt[j];
        			h_weighted_moving_sum[j] -= h_moving_sum[j] * h_min_weight;
        			h_moving_sum[j] -= h_moving_window.getFirst()[j];
        		}
        		h_moving_window.removeFirst();
        		if(0 <= i - h_r) {
        			double[] tmp = h_rslt.removeFirst();
        			for(int j = 0; j < src.length; j++) tmp[j] += h_moving_sum[j] * h_max_weight - h_weighted_moving_sum[j];
        			rslt.add(tmp);
        		}
			}
        }
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
	public static class Sobel {
		ArrayList<int[]> x_h;
		ArrayList<int[]> y_h;
		ArrayList<int[][]> rslt_list;
		ArrayList<int[][]> gradient_tmp;
		public void Add(int frame, int[] grayscale, boolean last) {
			int[] x_v = new int[grayscale.length];
			int[] y_v = new int[grayscale.length];
			x_v[0] = grayscale[1] - grayscale[0];
			y_v[0] = grayscale[1] + grayscale[0] * 3;
			for(int i = 1; i < grayscale.length - 1; i++) {
				x_v[i] = grayscale[i + 1] - grayscale[i - 1];
				y_v[i] = grayscale[i - 1] + grayscale[i] * 2 + grayscale[i + 1];
			}
			x_v[grayscale.length - 1] = grayscale[grayscale.length - 1] - grayscale[grayscale.length - 2];
			y_v[grayscale.length - 1] = grayscale[grayscale.length - 2] + grayscale[grayscale.length - 1] * 3;		
			if(frame == 0) {
				x_h = new ArrayList<>();
				y_h = new ArrayList<>();
				rslt_list = new ArrayList<>();
				gradient_tmp = new ArrayList<>();
				x_h.add(x_v);
				x_h.add(x_v);
				y_h.add(y_v);
				y_h.add(y_v);
			}
			
			for(int i = frame; i <= frame + (last ? 1 : 0); i++) {
				x_h.add(x_v);
				y_h.add(y_v);
				int[][] gradient = new int[grayscale.length][2];
				for(int j = 0; j < grayscale.length; j++) {
					int x_sum = x_h.get(0)[j] + x_h.get(1)[j] * 2 + x_h.get(2)[j];
					int y_sum = y_h.get(2)[j] - y_h.get(0)[j];
					if(x_sum == y_sum) gradient[j][0] = gradient[j][1] = 0;
					else {
						double rad = Math.atan2(y_sum, x_sum);
		        		if(rad < 0) rad += Math.PI * 2;
		        		gradient[j][0] = (int)Math.floor((rad + Math.PI / 4) / (Math.PI / 2)) % 4;
		        		gradient[j][1] = (int)Math.sqrt(x_sum * x_sum + y_sum * y_sum);
					}
				}
				x_h.remove(0);
				y_h.remove(0);
				rslt_list.add(gradient);
				
				/
				gradient_tmp.add(gradient);
				if(gradient_tmp.size() < 2) continue;
				int[][] canny = new int[grayscale.length][2];
				if(gradient_tmp.size() == 3) {
					for(int j = 1; j < grayscale.length - 1; j++) {
						
						if((gradient_tmp.get(1)[j][0] % 2 == 0 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(1)[j - 1][1], gradient_tmp.get(1)[j + 1][1])) ||
							(gradient_tmp.get(1)[j][0] % 2 == 1 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j][1], gradient_tmp.get(2)[j][1]))
	        			) {
							canny[j][0] = gradient_tmp.get(1)[j][0];
							canny[j][1] = Math.min(1, gradient_tmp.get(1)[j][1]);
						}
						
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

						/
						//int sector = (int)Math.floor((gradient_tmp.get(1)[j][0] + Math.PI / 8) / (Math.PI / 4)) % 8;
						if((sector % 4 == 0 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(1)[j - 1][1], gradient_tmp.get(1)[j + 1][1])) ||
							(sector % 4 == 1 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j - 1][1], gradient_tmp.get(2)[j + 1][1])) ||
							(sector % 4 == 2 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j][1], gradient_tmp.get(2)[j][1])) ||
							(sector % 4 == 3 && gradient_tmp.get(1)[j][1] >= Math.max(gradient_tmp.get(0)[j + 1][1], gradient_tmp.get(2)[j - 1][1]))
	        			) {
							canny[j][0] = gradient_tmp.get(1)[j][0];
							canny[j][1] = Math.min(1, gradient_tmp.get(1)[j][1]);
						}
						/


					}
					gradient_tmp.remove(0);
				}
				rslt_list.add(canny);
				/
				
				
			}
		}
	}
	*/
	
	public static class HarmonicSeperation {
		int radius;
		ArrayList<IndexValuePair>[] median;
		ArrayList<double[]> rslt;
		public HarmonicSeperation(int radius){
			this.radius = radius;
			rslt = new ArrayList<>();
		}
		public void Add(int frame, double[] value, boolean last) {
			if(frame == 0) {
				median = new ArrayList[value.length];
				for(int i = 0; i < value.length; i++) median[i] = new ArrayList<>();
			}
			for(int i = frame + (frame == 0 ? -radius : 0); i <= frame + (last ? radius : 0); i++) {
				for(int j = 0; j < value.length; j++) {
					int low = 0;
					int high = median[j].size();
					while(low < high) {
						int mid = (low + high) / 2;
						if(value[j] < (double)median[j].get(mid).value) high = mid;
						else low = mid + 1;
					}
					median[j].add(low, new IndexValuePair(i, value[j]));
				}
				if(i < radius) continue;
				double[] tmp = new double[value.length];
				for(int j = 0; j < value.length; j++) {
					tmp[j] = (double)median[j].get(radius).value; //(double)median[j].get(radius).value;
					for(int k = median[j].size() - 1; k >= 0; k--) {
						if(median[j].get(k).index == i - radius * 2) {
							median[j].remove(k);
							break;
						}
					}
				}
				rslt.add(tmp);
			}
        }
	}
	
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
		
		try {
			long stamp = System.currentTimeMillis();
			
			
			
			
			
			int keyStart = 0;
			int keyCount = 88;
			//double[] keyFrequency = new double[88];
			//for(int i = 0; i < 88; i++) noteFrequency[i] = Math.pow(2, i / 12.) * A0;
			int sampleRateOrigin = 44100;
			double sampleRate = sampleRateOrigin;
			double[] kaiserWindowedSinc = new double[] { 1 };
			int decimation = 1;
			while((sampleRateOrigin * 0.5 / (decimation + 1) - A0 * Math.pow(2, (keyStart + keyCount - 0.5) / 12)) >= sampleRateOrigin * 0.01) decimation++;
			//decimation = 1;
			if(decimation > 1) {
				sampleRate /= decimation;
				kaiserWindowedSinc = kaiserWindowedSinc(A0 * Math.pow(2, (keyStart + keyCount - 0.5) / 12) / sampleRateOrigin, sampleRate * 0.5 / sampleRateOrigin, 80);
			}
			System.out.println("kaiser N : " + kaiserWindowedSinc.length);	
			System.out.println(sampleRate);
			

            double frame_sec = 0.01;
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
            	
            	win[i] = 0.5 * (1 - Math.cos(2.0 * Math.PI * i / (win.length - 1)));
            	//win[i] = 1;
            	
            	win_sum += win[i];
            }
            for(int i = 0; i < win.length; i++) win[i] *= 1 / win_sum;


			/*
			double frame_sec = 0.02;
            int frame_len = (int)(sampleRate * frame_sec);
            double[] win = new double[frame_len * 10];
            double win_sum = 0;
    		{
    			int overlap = (win.length - frame_len) / 2;
    			for(int i = 0; i < win.length; i++) {
    				if(i < overlap) win[i] = Math.sin(Math.PI * i / (2 * (overlap - 1)));
    				else if(i < overlap + frame_len) win[i] = 1;
    				else win[i] = win[win.length - 1 - i];
    				
    				System.out.println(i + " : " + win[i]);
        			
        			win_sum += win[i];
    			}
    		}
    		for(int i = 0; i < win.length; i++) win[i] *= 2 / win_sum;
    		*/


    		
            
            int n_fft = 2;
            int mel_compression = 2;
            
            while(n_fft < win.length || A0 * (Math.pow(2, (keyStart + 0.5) / 12) - Math.pow(2, (keyStart - 0.5) / 12)) < sampleRate / n_fft) n_fft <<= 1;
			//while(n_fft < win.length || A0 * (Math.pow(2, (keyStart - 0.5 + 1. / mel_compression) / 12) - Math.pow(2, (keyStart - 0.5) / 12)) < sampleRate / n_fft) n_fft <<= 1;
            
            //n_fft = 8192;
            
            n_fft = n_fft * 2 * 2 * 2;
            System.out.println("1px : " + (sampleRate / n_fft) + "Hz");
            
            System.out.println("win : " + win.length + ", n_fft : " + n_fft);
            
            
            double[] sampleBuffer = new double[win.length];
            int sampleOffset = sampleBuffer.length - frame_len;
            
            int mel_px = (int)Math.floor(2 * n_fft / win_sum);
            int[][] keyRange = new int[keyCount][2];
            
            double equal_loudness_max = 0;
            

            double[] eq = new double[(int)(Math.min(A0 * Math.pow(2, (keyStart + keyCount - 0.5) / 12), sampleRate * 0.5) / sampleRate * n_fft)];
            System.out.println(eq.length);
            for(int i = 0; i < eq.length; i++) {
            	double pivot = 1000.;
            	double freq = (double)i * sampleRate / n_fft;
            	double h_pivot = ((1037918.48 - pivot * pivot) * (1037918.48 - pivot * pivot) + 1080768.16 * pivot * pivot) / ((9837328 - pivot * pivot) * (9837328 - pivot * pivot) + 11723776 * pivot * pivot);
            	double n_pivot = (pivot / (6.8966888496476 * Math.pow(10, -5))) * Math.sqrt(h_pivot / ((pivot * pivot + 79919.29) * (pivot * pivot + 1345600)));
            	double h_freq = ((1037918.48 - freq * freq) * (1037918.48 - freq * freq) + 1080768.16 * freq * freq) / ((9837328 - freq * freq) * (9837328 - freq * freq) + 11723776 * freq * freq);
            	double n_freq = (freq / (6.8966888496476 * Math.pow(10, -5))) * Math.sqrt(h_freq / ((freq * freq + 79919.29) * (freq * freq + 1345600)));
            	eq[i] = Math.abs(n_freq / n_pivot);
            	//System.out.println(i + " : " + eq[i]);
            	if(equal_loudness_max < eq[i]) equal_loudness_max = eq[i];
            }


            /*
            double[] hz = new double[] { 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500 };
            double[] af = new double[] { 0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349, 0.330, 0.315, 0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250, 0.246, 0.244, 0.243, 0.243, 0.243, 0.242, 0.242, 0.245, 0.254, 0.271, 0.301 };
            double[] lu = new double[] { -31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3, -8.1, -6.2, -4.5, -3.1, -2.0, -1.1, -0.4, 0.0, 0.3, 0.5, 0.0, -2.7, -4.1, -1.0, 1.7, 2.5, 1.2, -2.1, -7.1, -11.2, -10.7, -3.1 };
            double[] tf = new double[] { 78.5, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9, 14.4, 11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3, -4.2, -6.0, -5.4, -1.5, 6.0, 12.6, 13.9, 12.3 };
            double[] eq = new double[(int)(Math.min(A0 * Math.pow(2, (keyStart + keyCount - 0.5) / 12), hz[hz.length - 1]) / sampleRate * n_fft)];
            for(int i = 0; i < eq.length; i++) {
            	double freq = i * sampleRate / n_fft;
            	int low = 0;
            	int high = hz.length;
            	while(low < high) {
            		int mid = (low + high) / 2;
            		if(freq < hz[mid]) high = mid;
            		else low = mid + 1;
            	}
            	if(0 < low) {
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

            
            System.out.println("eq_max : " + equal_loudness_max);
            
            //for(int i = 0; i < eq.length; i++) eq[i] /= equal_loudness_max;
            
            
            /*
            for(int i = 0; i < keyCount; i++) {
            	System.out.println(i + " : " + keyRange[i][0] + ", " + keyRange[i][1]);
            }
            */

            /*
            BufferedImage img = new BufferedImage(eq.length, 500, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < img.getWidth(); i++) {
            	img.setRGB(i, (int)((img.getHeight() - 1) * (1 - eq[i])), 0xFFFFFFFF);
            }
            ImageIO.write(img, "PNG", new File("output_eq.png"));//
            */
            
            
            /*
            FileOutputStream smp = new FileOutputStream(new File("sample2.pcm"));
            double[] freqArray = new double[] { 1046.502, 523.2511, 880.0 };
            
            for(int i = 0; i < freqArray.length; i++) {
            	double sec = 1;
            	double period = sampleRate / freqArray[i];
            	
            	for(int j = 0; j < (int)(sec * sampleRate); j++) {
                	double angle = 2.0 * Math.PI * j / period;
                	short val = (short)(Math.sin(angle) * 32767  * 0.5);
                	smp.write(val & 0xFF);
                	smp.write((val >> 8) & 0xFF);
            	}
            	
            }
            FileInputStream fis = new FileInputStream(new File("sample2.pcm"));
            */
            

            FileInputStream fis = new FileInputStream(new File("sample/sample4.pcm"));
            boolean mel_gaus = true;
            boolean toDB = false;
            double dbRange = 80;
            boolean apply_eq = false;
            //double mel_radius = 0.75;
            double mel_radius = 0.5 / mel_compression;
            boolean is_zscore_threshold = false; //이게 true가 맞는가, false가 맞는가 확인해야한다.
            double zscore_threshold = 3;
            //adaptive
            double c = 5;
            
            int mel_size = (int)(12 * Math.log((eq.length - 1) * sampleRate / n_fft / 55) / Math.log(2) + 9 - mel_radius);
            
            double max_db = 0;

            
            /*
            for(int i = 0; i < 500; i++) {
            	int k = (int)Math.round(12 * mel_compression * Math.log(1 + 1.0 * i) / Math.log(2));
            	System.out.println(k);
            	//이전과 비교했을때, 1이하의 차이라면 탈출시켜야해
            }
            */
            
            /*
            for(int i = 0; i < mel_range.length; i++) {
            	double n = (double)i / mel_compression;
            	double center = Math.pow(2, (n - 9) / 12) * 55;
            	mel_range[i][0] = Math.max(0, (int)Math.floor(center * Math.pow(2, -mel_radius / 12) * n_fft / sampleRate));
            	mel_range[i][1] = Math.min(eq.length - 1, (int)Math.floor(center * Math.pow(2, mel_radius / 12) * n_fft / sampleRate));
            	
            	
            	//System.out.println(i + " : " + mel_range[i][0] + " ~ " + mel_range[i][1]);
            	
            	mel_weight[i] = new double[mel_range[i][1] - mel_range[i][0] + 1];
            	double weight_sum = 0;
            	for(int j = mel_range[i][0]; j <= mel_range[i][1]; j++) {
            		double note1 = 12 * Math.log((j + 0.) * sampleRate / n_fft / 55) / Math.log(2) + 9;
            		double note2 = 12 * Math.log((j + 1.) * sampleRate / n_fft / 55) / Math.log(2) + 9;
            		double weight = (mel_radius - (Math.max(n, note1) - Math.min(n, note2))) / mel_radius;
            		mel_weight[i][j - mel_range[i][0]] = weight;
            		weight_sum += weight;
            	}
            	for(int j = 0; j < mel_weight[i].length; j++) mel_weight[i][j] /= weight_sum;
            }
            */

            
            double bin_Hz = sampleRate / n_fft;
            double mel_Hz = 2 * sampleRate / win_sum; //이게 피아노 음 하나의 Hz공식임, 밑에꺼랑 합쳐서 계산해야된다.

            
            PriorityQueue<IndexValuePair> pq = new PriorityQueue<>((pair1, pair2)->{ return Double.valueOf(pair2.value).compareTo(pair1.value); });
            
            //System.out.println((int)Double.NEGATIVE_INFINITY);
            
            /*
            int[][] mel_range = new int[keyCount * mel_compression][2];
            for(int i = 0; i < eq.length; i++) {
            	int start = Math.max(0, (int)Math.floor((12 * Math.log(i * sampleRate / n_fft / A0) / Math.log(2) - keyStart + 0.5) * mel_compression));
            	int end = Math.min(mel_range.length - 1, Math.max(-1, (int)Math.floor((12 * Math.log((i + 1) * sampleRate / n_fft / A0) / Math.log(2) - keyStart + 0.5) * mel_compression)));
            	for(int j = start; j <= end; j++) {
            		if(mel_range[j][0] == 0) mel_range[j][0] = i;
            		mel_range[j][1] = i;
            	}
            }
            */
            int mel_comp = 1;
            
            double[] mel_lerp = new double[keyCount * mel_comp];
            for(int i = 0; i < mel_lerp.length; i++) mel_lerp[i] = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0 / bin_Hz;

            

            int[][] mel_range = new int[keyCount * mel_comp][2];
            double[] mel_freq = new double[mel_range.length];
            double[][] mel_weight = new double[mel_range.length][];
            
            for(int i = 0; i < mel_range.length; i++) {
            	
            	mel_freq[i] = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0;
            	
            	double center = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0 / bin_Hz;
            	
            	
            	//double center = Math.pow(2, ((i - keyStart - mel_comp / 2) / (double)mel_comp) / 12) * A0;
            	/*
            	double left = center * Math.pow(2, -1. / 12);
            	double right = center * Math.pow(2, 1. / 12);
            	*/
            	
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
            
            
            double gaus_max = 0;
            
            double prev_s = 0;
            double energ = 0;
            int dbCompression = 50;
            
            //double[][] note_range = new double[eq.length][2];
            //double note_range_sum = 0;
            //double note_range_max = 0;

            double max_energy = 0;
            
            ArrayList<Boolean> onset_list = new ArrayList<>();

            ArrayList<double[]> dbfs_list = new ArrayList<>();
            ArrayList<double[]> dbfs2_list = new ArrayList<>();
            
            ArrayList<int[]> grayscale_list = new ArrayList<>();
            ArrayList<int[]> morph_list = new ArrayList<>();
            ArrayList<int[]> tmp_list = new ArrayList<>();
            

            ArrayList<Double> smooth_energy_list = new ArrayList<>();
            /*
            ArrayList<double[]> sum_list = new ArrayList<>();
            double[] sum_total = new double[mel_size * mel_compression];
            ArrayList<double[]> sumsq_list = new ArrayList<>();
            double[] sumsq_total = new double[mel_size * mel_compression];
            ArrayList<double[][]> sumxy_list = new ArrayList<>();
            double[][] sumxy_total = new double[mel_size * mel_compression][mel_size * mel_compression];
            
            for(int i = 0; i < correlation_w_r * 2 + 1; i++) {
            	sum_list.add(new double[mel_size * mel_compression]);
            	sumsq_list.add(new double[mel_size * mel_compression]);
            	sumxy_list.add(new double[mel_size * mel_compression][mel_size * mel_compression]);
            }
            */
            
            ArrayList<int[][]> correlation_list = new ArrayList<>();
            
            ArrayList<int[]> sum_of_harmonic_list = new ArrayList<>();
            
            ArrayList<boolean[]> valid_list = new ArrayList<>();
            
            double max_zscore = 0;
            
            double prev_val = 0;
            
            //double[] prev = new double[mel_size * mel_compression];
            //double[] prev = new double[eq.length];
            
            
            double max_corr_sum = 0;
            ArrayList<int[]> corr_sum_list = new ArrayList<>();

            /*
            ArrayList<Double> filteredY_list = new ArrayList<>();
            ArrayList<Integer> signal_list = new ArrayList<>();
            double sum = 0;
            double sumsq = 0;
            for(int i = 0; i < 10; i++) {
            	filteredY_list.add(0.);
            }
            */

            //double[] sum = new double[prev.length];
            //double[] sumsq = new double[prev.length];
            ArrayList<double[]> note_list = new ArrayList<>();
            

            
            int[] prev_grayscale = new int[mel_size * mel_compression];
            
            ArrayList<int[]> freq_list = new ArrayList<>(); 
            
            ArrayList<double[][]> gradient_list = new ArrayList<>();

            
            long t_stamp = System.currentTimeMillis();
            
            //double[][] mel_range = new double[mel_size * mel_compression][2];
            
            
            double[] prev = new double[n_fft / 2];
            
            
            /*
            int[][] mag_range = new int[eq.length][2];
            double mag_radius = 0.5;
            
            for(int j = 0; j < eq.length; j++) {
            	double freq1 = (j + 0.) * sampleRate / n_fft;
            	double freq2 = (j + 1.) * sampleRate / n_fft;
            	
            	mag_range[j][0] = Math.max(0, j - mel_px / 2);
            	mag_range[j][1] = Math.min(eq.length - 1, j + mel_px / 2);
				//mag_range[j][0] = Math.max(0, (int)Math.floor(Math.min(freq1 - mel_Hz * mag_radius, freq1 * Math.pow(2, -mag_radius / 12)) * n_fft / sampleRate));
				//mag_range[j][1] = Math.min(eq.length - 1, (int)Math.floor(Math.max(freq2 + mel_Hz * mag_radius, freq2 * Math.pow(2, mag_radius / 12)) * n_fft / sampleRate));
            }
            */
            
            /*
            double dog_radius = 1.0;
            int[][] dog_range = new int[mel_lerp.length][2];
            for(int i = 0; i < mel_lerp.length; i++) {
            	double n = (double)i / mel_compression;
            	double center = Math.pow(2, (n - 9) / 12) * 55;
            	dog_range[i][0] = (int)Math.round((12 * Math.log((center - mel_Hz * dog_radius) / 55.0) / Math.log(2) + 9) * mel_compression);
            	dog_range[i][1] = (int)Math.round((12 * Math.log((center + mel_Hz * dog_radius) / 55.0) / Math.log(2) + 9) * mel_compression);
            	dog_range[i][0] = Math.max(0, Math.min(mel_lerp.length - 1, dog_range[i][0]));
            	dog_range[i][1] = Math.max(0, Math.min(mel_lerp.length - 1, dog_range[i][1]));
            }
            */

            Sobel sobel = new Sobel();
            //Prewitt prewitt = new Prewitt(7);
            GausFilter gausFilter1 = new GausFilter(mel_px / 2, mel_px / 2);
            GausFilter gausFilter2 = new GausFilter(mel_px * 3, mel_px * 3);
            
            HarmonicSeperation hs = new HarmonicSeperation(0);
            /*
            int start = -1;
            int peak = -1;
            ArrayList<Peak> peaks = new ArrayList<>();
			ArrayList<Integer> set = new ArrayList<>();
			*/
            ArrayList<Integer> max_idx_list = new ArrayList<>(); 
            
            ArrayList<Double> sum_list = new ArrayList<>();
            ArrayList<Double> sumsq_list = new ArrayList<>();
            int frames = 10;
            double frames_sum = 0;
            double frames_sumsq = 0;
            for(int i = 0; i < frames; i++) {
            	sum_list.add(0.);
            	sumsq_list.add(0.);
            }

            double comp_max = 0;
            
            double a = 90. / 180. * Math.PI; //160.
            
            double minDB = 20 * Math.log10(1. / 65536);
            
            
            int scale = 10;

            int frame = 0;//임시
            int frame2 = 0;

            //System.out.println((6.02 * 13 + 1.76) + "DB");
            
            /*
            int[][] mel_range = new int[keyCount * mel_compression][2];
            for(int i = 0; i < mel_range.length; i++) {
            	double center = Math.pow(2, (double)i / mel_compression / 12) * A0;
            	
            	//System.out.println(center);
            	
            	mel_range[i][0] = Math.max(0, (int)Math.floor(center * Math.pow(2, -mel_radius / 12) * n_fft / sampleRate));
            	mel_range[i][1] = Math.min(eq.length - 1, (int)Math.floor(center * Math.pow(2, mel_radius / 12) * n_fft / sampleRate));
            }
            */

            
            double[] prev_mel_lerp = new double[eq.length];
            
            int[] cdf = new int[256];
            int[] prev_cdf = cdf;
            int cdf_size = 4;
            
            int min_clip_limit = (int)Math.ceil((double)cdf_size * eq.length / cdf.length);
            double norm_clip_limit = 1;
            int clipLimit = min_clip_limit + (int)Math.round(norm_clip_limit * (cdf_size * eq.length - min_clip_limit));
            cdf_size = cdf_size <= 2 ? 2 : (cdf_size % 2 == 1 ? cdf_size - 1 : cdf_size);
            
            //System.out.println(cdf_size);
            List<double[]> list = new ArrayList<>();
            List<Double> max_list = new ArrayList<>();
            //int[] sampleBuffer = new int[(int)(44100 * 0.05) / 2 * 2];
            //int[] prevBuffer = new int[sampleBuffer.length / 2];
            
            
            int channelCount = 1;
            int div = channelCount * 2;
            byte[] buffer = new byte[4096];
            int bufferSize = 0;
            
            double mag_max = 0;

            
            int mel_comp2 = 1;
            int[][] mel_range2 = new int[keyCount * mel_comp2][2];
            for(int i = 0; i < n_fft / 2; i++) {
            	int key = (int)Math.round(Math.max(-1, 12 * Math.log((double)i * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp2);
            	System.out.println(Math.max(-1, 12 * Math.log((double)i * sampleRate / n_fft / A0) / Math.log(2)) * mel_comp2);
            	if(0 <= key && key < mel_range2.length) {
            		if(mel_range2[key][0] == 0) mel_range2[key][0] = i;
            		mel_range2[key][1] = i;
            	}
            }
            

            int[][] mag_range = new int[eq.length][2];
            double mag_radius = 0.5;
            for(int j = 0; j < eq.length; j++) {
            	double freq1 = (j + 0.) * sampleRate / n_fft;
            	double freq2 = (j + 1.) * sampleRate / n_fft;
            	
            	//mag_range[j][0] = Math.max(0, j - mel_px / 2);
            	//mag_range[j][1] = Math.min(eq.length - 1, j + mel_px / 2);
				//mag_range[j][0] = Math.max(0, (int)Math.floor(Math.min(freq1 - mel_Hz * mag_radius, freq1 * Math.pow(2, -mag_radius / 12)) * n_fft / sampleRate));
				//mag_range[j][1] = Math.min(eq.length - 1, (int)Math.floor(Math.max(freq2 + mel_Hz * mag_radius, freq2 * Math.pow(2, mag_radius / 12)) * n_fft / sampleRate));
            	mag_range[j][0] = Math.max(0, (int)Math.floor(freq1 * Math.pow(2, -mag_radius / 12) * n_fft / sampleRate));
            	mag_range[j][1] = Math.min(eq.length - 1, (int)Math.floor(freq2 * Math.pow(2, mag_radius / 12) * n_fft / sampleRate));
            }

            
            ArrayList<ArrayList<Peak>> peaks_list = new ArrayList<>();
            
            FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            
            List<Double> frame_list = new ArrayList<>();
            
            List<List<Mel>> mel_frames = new ArrayList<>();
            
            List<Integer> samp_list = new ArrayList<>(kaiserWindowedSinc.length);
            for(int i = 0; i < kaiserWindowedSinc.length / 2; i++) samp_list.add(0);
            int curr = decimation - 1;
            while((bufferSize = fis.read(buffer, 0, buffer.length)) != -1) {
            	for(int i = 0; i < bufferSize / div; i++) {
            		int samp = 0;
            		for(int j = 0; j < channelCount; j++) samp += (buffer[i * div + j * 2 + 1] << 8) | (buffer[i * div + j * 2] & 0xFF);
            		samp_list.add(samp / channelCount);
            		if(samp_list.size() < kaiserWindowedSinc.length) continue;
            		if(++curr >= decimation) {
            			curr = 0;
            			double sum = 0;
            			for(int j = 0; j < kaiserWindowedSinc.length; j++) sum += samp_list.get(j) * kaiserWindowedSinc[j];
            			//sampleBuffer[sampleOffset++] = sum;
            			sampleBuffer[sampleOffset++] = sum / 32767.;
            			
                        if(sampleOffset == sampleBuffer.length) {
                        	Complex[] complex = new Complex[n_fft];
        					for(int j = 0; j < n_fft; j++) complex[j] = new Complex(j < sampleBuffer.length ? sampleBuffer[j] * win[j] : 0, 0);
        					complex = fft.transform(complex, TransformType.FORWARD);

        					double[] mag = new double[n_fft / 2];
        					double[] phase = new double[n_fft / 2];
        					for(int j = 0; j < mag.length; j++) {
        						mag[j] = complex[j].abs(); //이건가 아래껀가 고민됨
        						//mag[j] = complex[j].abs();
        						
        						phase[j] = Math.atan2(complex[j].getImaginary(), complex[j].getReal());
        						
        						if(mag[j] < mag_max) mag_max = mag[j];
        					}
        					
        					/*
        					double[] mel = new double[mel_range2.length];
        					int[] freq = new int[mel_range2.length];
        					double[] phase = new double[mel_range2.length];
        					for(int j = 0; j < mel_range2.length; j++) {
        						for(int k = mel_range2[j][0]; k <= mel_range2[j][1]; k++) {
        							if(mel[j] < complex[k].abs()) {
        								freq[j] = k;
        								mel[j] = complex[k].abs();
        								phase[j] = Math.atan2(complex[k].getImaginary(), complex[k].getReal());
        							}
        						}
        						if(mag_max < mel[j]) mag_max = mel[j];
        					}
        					
        					
        					*/
        					

        					int[] curr_range = new int[] {0, 0};
    			            pq.add(new IndexValuePair(0, mag[0]));
    						for(int j = 0; j < eq.length; j++) {
    							//System.out.println(mag_range[j][0] + " ~ " + mag_range[j][1]);
    							for(int k = curr_range[1] + 1; k <= mag_range[j][1]; k++) pq.add(new IndexValuePair(k, mag[k]));
    							curr_range = mag_range[j];
    							while(pq.peek().index < curr_range[0]) pq.remove();
    							mag[j] = pq.peek().value;
    							//dbfs[j] = Math.max(0, 20 * Math.log10(pq.peek().value) + dbNorm) / dbNorm;
    							//grayscale[j] = (int)(Math.max(0, 20 * Math.log10(dbfs[j]) + dbNorm) / dbNorm * 0xFF);
    						}
    						pq.clear();
        					
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
									peaks.add(new Peak(frame, start, peak, end, value[peak[0]]));
    								//sum += value;
    								//if(comp_max < value) comp_max = value;
    							}
    						}
							peaks_list.add(peaks);
							
							
							dbfs_list.add(mag);
        					dbfs2_list.add(phase);
        					
        					//mag = blur(mag, 50);
        					
        					/*
        					double[] mel = new double[mel_range.length];
        					for(int j = 0; j < mel.length; j++) {
        						for(int k = mel_range[j][0]; k <= mel_range[j][1]; k++) {
        							mel[j] += mag[k] * mel_weight[j][k - mel_range[j][0]];
        						}
        					}
        					
        					mag = mel;
        					*/
        					
        					/*
        					for(int j = 0; j < mag.length; j++) {
        						if(mag_max < mag[j]) mag_max = mag[j];
        					}
        					dbfs_list.add(mag);
        					dbfs2_list.add(angle);
        					*/


        					sampleOffset = sampleBuffer.length - frame_len;
                            System.arraycopy(sampleBuffer, frame_len, sampleBuffer, 0, sampleOffset);
                            frame++;
                        }
            		}
            		samp_list.remove(0);
            	}
            }
            System.out.println(System.currentTimeMillis() - stamp);
            
            

            for(int i = 0; i < peaks_list.size(); i++) {
            	List<Peak> peaks = peaks_list.get(i);
            	
            	List<Mel> mels = new ArrayList<>();
            	double[] mag = dbfs_list.get(i);
            	double[] angle = dbfs2_list.get(i);
            	
            	for(int j = 0; j < peaks.size(); j++) {
            		int peak_index = (peaks.get(j).peak[0] + peaks.get(j).peak[1]) / 2;
            		mels.add(new Mel(peak_index * bin_Hz, mag[peak_index] / mag_max, angle[peak_index]));
            	}
            	mel_frames.add(mels);
            }

            
            /*
            BufferedImage dbfs = new BufferedImage(dbfs_list.size(), dbfs_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < dbfs_list.size(); i++) {
            	double[] mel = dbfs_list.get(i);
            	double[] phase = dbfs2_list.get(i);
            	int[] freq = freq_list.get(i);
            	
            	List<Mel> mels = new ArrayList<>();
            	
            	for(int j = 1; j < mel.length - 1; j++) {

            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) mels.add(new Mel(mel_freq[j], mag[j] / mag_max));
            		//mag[j] > Math.max(mag[j - 1], mag[j + 1])
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1]) && j * bin_Hz > 20)
            		if(mel[j] > Math.max(mel[j - 1], mel[j + 1])) mels.add(new Mel(freq[j] * bin_Hz, mel[j] / mag_max, phase[j]));
            		
            		//mels.add(new Mel(Math.pow(2, j / 12.) * A0, mel[j] / mag_max, phase[j]));
            		
            		int col = (int)(mel[j] / mag_max * 0xFF);
            		dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(dbfs_list.get(i)[j] > Math.max(dbfs_list.get(i)[j - 1], dbfs_list.get(i)[j + 1])) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            	mel_frames.add(mels);
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            */
            
            byte[] pcmBytesArray = mel_to_pcm((int)sampleRate, frame_sec, mel_frames);
            
            pcmToWav((int)sampleRate, channelCount, 16, pcmBytesArray, "output.wav");
            
            
            /*
            BufferedImage dbfs = new BufferedImage(dbfs_list.size(), dbfs_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < dbfs_list.size(); i++) {
            	double[] mag = dbfs_list.get(i);
            	//double[] angle = dbfs2_list.get(i);
            	
            	List<Mel> mels = new ArrayList<>();
            	
            	for(int j = 1; j < mag.length - 1; j++) {

            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1])) mels.add(new Mel(mel_freq[j], mag[j] / mag_max));
            		//mag[j] > Math.max(mag[j - 1], mag[j + 1])
            		//if(mag[j] > Math.max(mag[j - 1], mag[j + 1]) && j * bin_Hz > 20) mels.add(new Mel(j * bin_Hz, mag[j] / mag_max, angle[j]));
            		
            		int col = (int)(mag[j] / mag_max * 0xFF);
            		dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            		//if(dbfs_list.get(i)[j] > Math.max(dbfs_list.get(i)[j - 1], dbfs_list.get(i)[j + 1])) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            	//mel_frames.add(mels);
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            */
            
            /*
            byte[] pcmBytesArray = mel_to_pcm((int)sampleRate, frame_sec, mel_frames);
            
            pcmToWav((int)sampleRate, channelCount, 16, pcmBytesArray, "output.wav");
            */
            
            /*
            int numColors = dbfs_list.get(0).length;
            int[] colors = new int[numColors];
            for (int i = 0; i < numColors; i++) {
            	//float hue = (float)i / (numColors - 1);
            	float hue = (float)i / (numColors - 1) * 0.8f;
                //float hue = (float)(i % 12) / 11 * 0.8f;
                colors[i] = Color.getHSBColor(hue, 1.0f, 1.0f).getRGB(); // 색상 생성
            }
            BufferedImage col_test = new BufferedImage(numColors, 1, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < col_test.getWidth(); i++) {
            	for(int j = 0; j < col_test.getHeight(); j++) {
            		col_test.setRGB(i, j, colors[i]);
            	}
            }
            ImageIO.write(col_test, "PNG", new File("output_col_test.png"));//
            */

            /*
            BufferedImage grayscale = new BufferedImage(grayscale_list.size(), grayscale_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < grayscale.getWidth(); i++) {
            	for(int j = 1; j < grayscale.getHeight() - 1; j++) {
            		int col = grayscale_list.get(i)[j];
        			

        			//if(dbfs_list.get(i)[j] > 0) col = 0xFF;
        			//else col = 0;
        			
        			if(col > 0) col = 0xFF;

            		grayscale.setRGB(i, grayscale_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
        			
        			//if(j == max_idx_list.get(i)) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
        			if(i < valid_list.size() && valid_list.get(i)[j]) grayscale.setRGB(i, grayscale_list.get(0).length - 1 - j, 0xFFFF0000);
        			//dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, colors[signal_list.get(i)[j]]);
        			//if(dbfs_list.get(i)[j] > Math.max(dbfs_list.get(i)[j - 1], dbfs_list.get(i)[j + 1]) && (int)(dbfs_list.get(i)[j] * 0xFF) > 0) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
        			//if(grayscale_list.get(i)[j] > Math.max(grayscale_list.get(i)[j - 1], grayscale_list.get(i)[j + 1])) grayscale.setRGB(i, grayscale_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(grayscale, "PNG", new File("output_grayscale.png"));//
            */
            
            
            /*
            BufferedImage dbfs = new BufferedImage(dbfs_list.size(), dbfs_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < dbfs.getWidth(); i++) {
            	for(int j = 1; j < dbfs.getHeight() - 1; j++) {
        			//int col = (int)(Math.min(1, dbfs_list.get(i)[j]) * 0xFF);
            		
            		
            		int col = (int)(Math.min(1, dbfs_list.get(i)[j] / mag_max) * 0xFF);
            		
            		//int col = (int)(dbfs_list.get(i)[j]);
            		//int col = (int)(Math.min(1, dbfs_list.get(i)[j]) * 0xFF);
        			

        			//if(col > 0) col = 0xFF;
        			//else col = 0;
        			
        			//if(col > 0) col = 0xFF;

        			dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
        			
        			//if(j == max_idx_list.get(i)) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
        			//if(i < valid_list.size() && valid_list.get(i)[j]) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
        			//dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, colors[signal_list.get(i)[j]]);
        			//if(dbfs_list.get(i)[j] > Math.max(dbfs_list.get(i)[j - 1], dbfs_list.get(i)[j + 1]) && (int)(dbfs_list.get(i)[j] * 0xFF) > 0) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
        			//if(dbfs_list.get(i)[j] > Math.max(dbfs_list.get(i)[j - 1], dbfs_list.get(i)[j + 1])) dbfs.setRGB(i, dbfs_list.get(0).length - 1 - j, 0xFFFF0000);
            	}
            }
            ImageIO.write(dbfs, "PNG", new File("output_dbfs.png"));//
            */
            

            /*
            BufferedImage energy = new BufferedImage(smooth_energy_list.size(), 500, BufferedImage.TYPE_INT_RGB);
            for(int i = 1; i < smooth_energy_list.size() - 1; i++) {
            	int height = (int)((energy.getHeight() - 1) * (Math.abs(smooth_energy_list.get(i)) / max_energy));
            	
            	//boolean peak = smooth_energy_list.get(i) > Math.max(smooth_energy_list.get(i - 1), smooth_energy_list.get(i + 1));
            	
            	for(int j = 0; j <= height; j++) {
            		//energy.setRGB(i, (energy.getHeight() - 1) - j, peak ? 0xFFFF0000 : 0xFFFFFFFF);
            		energy.setRGB(i, (energy.getHeight() - 1) - j, smooth_energy_list.get(i) < 0 ? 0xFFFF0000 : 0xFFFFFFFF);
            	}
            }
            ImageIO.write(energy, "PNG", new File("output_energy_list.png"));
            */
            
            /*
            BufferedImage ang = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            BufferedImage deg = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            
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
            BufferedImage peak = new BufferedImage(peaks_list.size(), dbfs_list.get(0).length, BufferedImage.TYPE_INT_RGB);
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
            			peak.setRGB(i, dbfs_list.get(0).length - 1 - k, col);
            		}
            	}
            }
            ImageIO.write(peak, "PNG", new File("output_peak.png"));//
            */
            
            /*
            BufferedImage ang = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            BufferedImage deg = new BufferedImage(gradient_list.size(), gradient_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            
            for(int i = 0; i < ang.getWidth(); i++) {
            	for(int j = 1; j < gradient_list.get(0).length - 1; j++) {
            		//int degree = 0xFF;
            		int degree = gradient_list.get(i)[j][1];
            		//if(gradient_list.get(i)[j][1] > 1. / 65536) degree = 0xFF;
            		
            		degree = 0xFF;
            		
        			int col = 0;
        			
        			//int sec = (int)Math.floor((gradient_list.get(i)[j][0] + Math.PI / 4) / (Math.PI / 2)) % 4;
        			
        			//double angle = gradient_list.get(i)[j][0] + a / 2;
        			//int sec = (int)(Math.floor(angle / Math.PI) * 2 + (Math.floor(angle % Math.PI / a) > 0 ? 1 : 0)) % 4;
        			int sec = gradient_list.get(i)[j][0];
        			//degree = 0xFF;
        			if(sec == 0) col = (degree << 24 | degree << 16);
            		else if(sec == 1) col = (degree << 24 | degree << 8);
            		else if(sec == 2) col = (degree << 24 | degree);
            		else if(sec == 3) col = (degree << 24 | degree << 16 | degree);
        			
        			
        			//if(gradient_list.get(i)[j][1] > 10) ang.setRGB(i, gradient_list.get(0).length - 1 - j, col);
        			if(gradient_list.get(i)[j][1] > 0) ang.setRGB(i, gradient_list.get(0).length - 1 - j, col);
        			
        			//if(i < valid_list.size() && valid_list.get(i)[j]) ang.setRGB(i, gradient_list.get(0).length - 1 - j, col);
        			
        			int d = Math.min(0xFF, gradient_list.get(i)[j][1]);
        			
        			//if(gradient_list.get(i)[j][1] > 0) System.out.println(gradient_list.get(i)[j][1]);
        			//if(d > 10) d = 0xFF;
        			deg.setRGB(i, gradient_list.get(0).length - 1 - j, (d << 24 | d << 16 | d << 8 | d));
            	}
            }
            ImageIO.write(ang, "PNG", new File("output_gradient_angle.png"));//
            ImageIO.write(deg, "PNG", new File("output_gradient_magnitude.png"));//
            */

            //int[][] mel = new int[dbfs_list.size()][keyCount];
            
            
            /*
            BufferedImage peak = new BufferedImage(peaks_list.size(), dbfs_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < peaks_list.size(); i++) {
            	ArrayList<Peak> peaks = peaks_list.get(i);
            	for(int j = 0; j < peaks.size(); j++) {
            		for(int k = peaks.get(j).start; k <= peaks.get(j).end; k++) {
            			int col = 0;
            			
            			//int degree = (int)(peaks.get(j).value / (scale - 1.) * 0xFF);
            			int degree = (int)(peaks.get(j).value / mag_max * 0xFF);
            			
            			if(peaks.get(j).peak[0] <= k && k <= peaks.get(j).peak[1]) col = 0xFFFF0000;
            			else if(k < peaks.get(j).peak[0]) col = (degree << 24 | degree << 8) ; //0xFF00FF00;
            			else if(k > peaks.get(j).peak[1]) col = (degree << 24 | degree);
            			peak.setRGB(i, dbfs_list.get(0).length - 1 - k, col);
            			
            			mel[i][peaks.get(j).peak[0] / mel_comp] = Math.max(mel[i][peaks.get(j).peak[0] / mel_comp], (int)(degree / 255.0 * 127));
            			
            			if(mel[i][peaks.get(j).peak[0] / mel_comp] > 50) System.out.println(mel[i][peaks.get(j).peak[0] / mel_comp]);
            		}
            	}
            }
            ImageIO.write(peak, "PNG", new File("output_peak.png"));//
            */


            /*
            BufferedImage ang2 = new BufferedImage(gradient2_list.size(), gradient2_list.get(0).length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < ang2.getWidth(); i++) {
            	for(int j = 0; j < gradient2_list.get(0).length; j++) {
            		int degree = gradient2_list.get(i)[j][1];
            		//if(degree > 0) degree = 0xFF;
        			int col = 0;
        			int sec = gradient2_list.get(i)[j][0];
        			
        			if(sec == 0) col = (degree << 24 | degree << 16);
            		else if(sec == 1) col = (degree << 24 | degree << 8);
            		else if(sec == 2) col = (degree << 24 | degree);
            		else if(sec == 3) col = (degree << 24 | degree << 16 | degree);
        			
        			if(gradient2_list.get(i)[j][1] > 0) ang2.setRGB(i, gradient2_list.get(0).length - 1 - j, col);
            	}
            }
            ImageIO.write(ang2, "PNG", new File("output_gradient_angle2.png"));//
            */
           
            
            /*
            for(int i = 0; i < dbfs_list.get(0).length; i++) {
            	BufferedImage note_img = new BufferedImage(dbfs_list.size(), 500, BufferedImage.TYPE_INT_RGB);
            	for(int j = 0; j < dbfs_list.size(); j++) {
            		
            		for(int l = note_img.getHeight() - 1; l >= (note_img.getHeight() - 1) * Math.max(0, 1 - dbfs_list.get(j)[i]); l--) {
            			note_img.setRGB(j, l, 0xFFFFFFFF);
        			}
            	}
            	ImageIO.write(note_img, "PNG", new File("fr/output_img" + i + ".png"));//
            }
            */
            
            /*
        	BufferedImage mel_img = new BufferedImage(mel.length, mel[0].length, BufferedImage.TYPE_INT_RGB);
            for(int i = 0; i < mel.length; i++) {
            	for(int j = 0; j < mel[0].length; j++) {
            		for(int k = keyRange[j][0]; k <= keyRange[j][1]; k++) {
            			mel[i][j] = Math.max(mel[i][j], (int)(grayscale_list.get(i)[k] / 255.0 * 127));
            		}
            		int col = mel[i][j];
        			mel_img.setRGB(i, mel[0].length - 1 - j, (col << 24 | col << 16 | col << 8 | col));
            	}
            }
            ImageIO.write(mel_img, "PNG", new File("output_mel.png"));//
            */


            /*
            FileOutputStream fos = new FileOutputStream(new File("output.midi"));
            int ppqn = 1;
            fos.write(new byte[] {
        		0x4D, 0x54, 0x68, 0x64,
        		0x00, 0x00, 0x00, 0x06,
        		0x00, 0x00, 0x00, 0x01,
        		(byte)(ppqn >> 8 & 0x7F), (byte)(ppqn & 0xFF),
        		0x4D, 0x54, 0x72, 0x6B
            });
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            //int ppqn_us = 100000;
            int ppqn_us = (int)Math.round(1000000.0 * frame_sec);
            int inst = 0; // 피아노 : 0, 바이올린 : 40, 글로켄슈필 : 9
            bos.write(new byte[] {
        		0x00, (byte)0xFF, 0x51, 0x03,
        		(byte)(ppqn_us >> 16 & 0xFF),
        		(byte)(ppqn_us >> 8 & 0xFF),
                (byte)(ppqn_us & 0xFF),
                0x00, (byte)0xC0, (byte)inst, 
            });
            int delta = 0;
            int[] note_state = new int[keyCount];
            for(int i = 0; i < mel.length; i++) {
            	for(int j = 0; j < mel[i].length; j++){
            		if(mel[i][j] != note_state[j]) {
            			boolean pass = true;
            			for (int k = 3; k >= 0; k--)
        		        {
        		            int value = delta >> (7 * k) & 0x7F | (k != 0 ? 0x80 : 0x00);
        		            if (value != 0x80) pass = false;
        		            if (!pass) bos.write(value);
        		        }
            			bos.write(0x90);
        				bos.write(j + 21 + keyStart);
        				//System.out.println(mag_list.get(i)[j] / mag_max);
        				bos.write(mel[i][j]);
            			delta = 0;
            			note_state[j] = mel[i][j];
            		}
            	}
            	delta++;
            }
            bos.write(new byte[] {
            		0x00, (byte)0xFF, 0x2F, 0x00
    		});
            byte[] track_chunk = bos.toByteArray();
            fos.write(new byte[] {
        		(byte)(track_chunk.length >> 24 & 0xFF),
        		(byte)(track_chunk.length >> 16 & 0xFF),
        		(byte)(track_chunk.length >> 8 & 0xFF),
        		(byte)(track_chunk.length & 0xFF)	
            });
            fos.write(track_chunk);
            fos.flush();
            fos.close();
            */
            int[][] mel = new int[frame_list.size()][1];
            
            for(int i = 0; i < frame_list.size(); i++) {
            	mel[i][0] = (int)(frame_list.get(i) / mag_max * 127);
            	System.out.println(i + " : " + mel[i][0]);
            }
            
            FileOutputStream fos = new FileOutputStream(new File("output.midi"));
            int ppqn = 1;
            fos.write(new byte[] {
        		0x4D, 0x54, 0x68, 0x64,
        		0x00, 0x00, 0x00, 0x06,
        		0x00, 0x00, 0x00, 0x01,
        		(byte)(ppqn >> 8 & 0x7F), (byte)(ppqn & 0xFF),
        		0x4D, 0x54, 0x72, 0x6B
            });
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            //int ppqn_us = 100000;
            int ppqn_us = (int)Math.round(1000000.0 * frame_sec);
            int inst = 54; // 피아노 : 0, 바이올린 : 40, 글로켄슈필 : 9
            bos.write(new byte[] {
        		0x00, (byte)0xFF, 0x51, 0x03,
        		(byte)(ppqn_us >> 16 & 0xFF),
        		(byte)(ppqn_us >> 8 & 0xFF),
                (byte)(ppqn_us & 0xFF),
                0x00, (byte)0xC0, (byte)inst, 
            });
            int delta = 0;
            int[] note_state = new int[keyCount];
            for(int i = 0; i < mel.length; i++) {
            	for(int j = 0; j < mel[i].length; j++){
            		if(mel[i][j] != note_state[j]) {
            			boolean pass = true;
            			for (int k = 3; k >= 0; k--)
        		        {
        		            int value = delta >> (7 * k) & 0x7F | (k != 0 ? 0x80 : 0x00);
        		            if (value != 0x80) pass = false;
        		            if (!pass) bos.write(value);
        		        }
            			bos.write(0x90);
        				//bos.write(j + 21 + keyStart);
            			bos.write(60);
        				//System.out.println(mag_list.get(i)[j] / mag_max);
        				bos.write(mel[i][j]);
            			delta = 0;
            			note_state[j] = mel[i][j];
            		}
            	}
            	delta++;
            }
            bos.write(new byte[] {
            		0x00, (byte)0xFF, 0x2F, 0x00
    		});
            byte[] track_chunk = bos.toByteArray();
            fos.write(new byte[] {
        		(byte)(track_chunk.length >> 24 & 0xFF),
        		(byte)(track_chunk.length >> 16 & 0xFF),
        		(byte)(track_chunk.length >> 8 & 0xFF),
        		(byte)(track_chunk.length & 0xFF)	
            });
            fos.write(track_chunk);
            fos.flush();
            fos.close();
            
            
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
}
