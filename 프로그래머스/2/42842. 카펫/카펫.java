class Solution {
    public int[] solution(int brown, int yellow) {
        int total = brown + yellow;
        for(int i = 3; i <= (int)Math.sqrt(total); i++){
            if(total % i != 0) continue;
            int w = Math.max(i, total / i);
            int h = Math.min(i, total / i);
            if(brown == 2 * w + 2 * h - 4) return new int[] { w, h };
        }
        return null;
    }
}