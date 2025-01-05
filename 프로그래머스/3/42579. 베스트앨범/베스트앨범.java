import java.util.*;
class Music{
    int number;
    String genre;
    int play;
    Music(int number, String genre, int play){
        this.number = number;
        this.genre = genre;
        this.play = play;
    }
}
class Solution {
    public int[] solution(String[] genres, int[] plays) {
        HashMap<String, Integer> map = new HashMap<>();
        List<Music> musicList = new ArrayList<>();
        for(int i = 0; i < genres.length; i++){
            map.put(genres[i], map.getOrDefault(genres[i], 0) + plays[i]);
            musicList.add(new Music(i, genres[i], plays[i]));
        }
        Collections.sort(musicList, (music1, music2)->{
            int play1 = map.get(music1.genre);
            int play2 = map.get(music2.genre);
            if(play1 != play2) return play2 - play1;
            if(music1.play != music2.play) return music2.play - music1.play;
            return music1.number - music2.number;
        });
        List<Integer> musicNumberList = new ArrayList<>();
        HashMap<String, Integer> map2 = new HashMap<>();
        for(int i = 0; i < musicList.size(); i++){
            String genre = musicList.get(i).genre;
            if(map2.getOrDefault(genre, 0) < 2){
                musicNumberList.add(musicList.get(i).number);
                map2.put(genre, map2.getOrDefault(genre, 0) + 1);
            }
        }
        int[] answer = new int[musicNumberList.size()];
        for(int i = 0; i < answer.length; i++) answer[i] = musicNumberList.get(i);
        return answer;
    }
}