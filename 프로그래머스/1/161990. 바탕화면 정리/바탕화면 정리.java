class Solution {
    public int[] solution(String[] wallpaper) {
        int[] answer = new int[4];
        answer[0] = wallpaper.length;
        answer[1] = wallpaper[0].length();
        answer[2] = 0;
        answer[3] = 0;
        for(int i = 0; i < wallpaper.length; i++){
            int min_y = wallpaper[i].indexOf("#");
            if(0 <= min_y){
                int max_y = wallpaper[i].lastIndexOf("#") + 1;
                answer[0] = Math.min(answer[0], i);
                answer[1] = Math.min(answer[1], min_y);
                answer[2] = Math.max(answer[2], i + 1);
                answer[3] = Math.max(answer[3], max_y);
            }
        }
        return answer;
    }
}