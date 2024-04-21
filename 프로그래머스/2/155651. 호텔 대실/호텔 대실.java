class Solution {
    public int solution(String[][] book_time) {
        int[] time = new int[24 * 60];
        for(int i = 0; i < book_time.length; i++){
            int start_hh = Integer.parseInt(book_time[i][0].substring(0, 2));
            int start_mm = Integer.parseInt(book_time[i][0].substring(3, 5));
            int end_hh = Integer.parseInt(book_time[i][1].substring(0, 2));
            int end_mm = Integer.parseInt(book_time[i][1].substring(3, 5)) + 9;
            if(end_mm >= 60){
                end_hh += end_mm / 60;
                end_mm %= 60;
            }
            for(int j = start_hh * 60 + start_mm; j <= Math.min(time.length - 1, end_hh * 60 + end_mm); j++) time[j]++;
        }
        int max = 0;
        for(int i = 0; i < time.length; i++) max = Math.max(max, time[i]);
        return max;
    }
}