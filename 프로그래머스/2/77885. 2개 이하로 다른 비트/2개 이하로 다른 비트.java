class Solution {
    public long[] solution(long[] numbers) {
        long[] answer = new long[numbers.length];
        for(int i = 0; i < numbers.length; i++){
            long tmp = numbers[i];
            int pos = 0;
            while((tmp & 1) != 0){
                tmp >>= 1;
                pos++;
            }
            answer[i] = numbers[i] | (1L << pos);
            if(pos > 0) answer[i] &= ~(1L << (pos - 1));
        }
        return answer;
    }
}