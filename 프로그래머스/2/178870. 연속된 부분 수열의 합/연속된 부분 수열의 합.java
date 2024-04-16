class Solution {
    public int[] solution(int[] sequence, int k) {
        int[] range = new int[2];
        range[0] = range[1] = sequence.length - 1;
        int sum = sequence[range[0]];
        while(sum != k){
            if(sum < k) sum += sequence[--range[0]];
            else sum -= sequence[range[1]--];
        }
        while(0 < range[0] && sequence[range[0] - 1] == sequence[range[1]]){
            range[0]--;
            range[1]--;
        }
        return range;
    }
}