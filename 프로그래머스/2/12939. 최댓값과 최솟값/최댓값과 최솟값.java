class Solution {
    public String solution(String s) {
        String[] sp = s.split("\\s");
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for(int i = 0; i < sp.length; i++){
            min = Math.min(min, Integer.parseInt(sp[i]));
            max = Math.max(max, Integer.parseInt(sp[i]));
        }
        return min + " " + max;
    }
}