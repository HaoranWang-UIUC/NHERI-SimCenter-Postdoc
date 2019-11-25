//  Question_1
//  Created by Haoran Wang on 11/23/19.
//

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>

using namespace std;
typedef map<string, int> StrIntMap;

void CountFiles(istream& in, unsigned long* line_count, unsigned long* word_count, unsigned long* uniq_word_count, unsigned long* palin_word_count, unsigned long* char_count, int i)
{

    string word;
    string line;
    StrIntMap word_table;
    line_count[i] = 0;
    word_count[i] = 0;
    uniq_word_count[i] = 0;
    palin_word_count[i] = 0;
    char_count[i] = 0;

    //Read each line and do counting
    while(getline(in,line))
    {
        size_t count_temp = line.length() - count(line.begin(), line.end(), ' ');//count non-space characters on each line
        
        if (count_temp > 0)//it's not an empty line
        {
            line_count[i]++;
            char_count[i] += count_temp;
            
            //Replace all the punctuations with space for word parsing
            //The 'replace' below should be modified based on the definition of 'word'
            replace(begin(line),end(line),',',' ');
            replace(begin(line),end(line),'.',' ');
            replace(begin(line),end(line),':',' ');
            
            istringstream iss(line);
            
            while (iss >> word)
            {
                word_count[i]++;
                ++word_table[word];//put words into hash table, update count
                if ( (word == string(word.rbegin(), word.rend())) && word_table[word]==1 )
                    palin_word_count[i]++;
            }
        }
    }
    
    uniq_word_count[i] += word_table.size();
    
}//end of CountFiles function
//******************************************************
//******************************************************
int main(int argc, char** argv) {

    if (argc < 3)//check command
        {cout << "Error: At least 1 input and 1 output are required!" << endl;
            return(EXIT_FAILURE);
        }
    
    unsigned long line_count[argc-2], word_count[argc-2], uniq_word_count[argc-2], palin_word_count[argc-2], char_count[argc-2];
    
    for (int i=0; i<argc-2; i++) //loop over each input file
    {
        ifstream input_file(argv[i+1]);
        
        if (!input_file)//check input files
        {cout << "Error: " << argv[i+1] << " could not be found!" << endl;
         exit(EXIT_FAILURE);
        }
        
        CountFiles(input_file, line_count, word_count, uniq_word_count, palin_word_count, char_count, i);
    }
    
    //write results to output file
    ofstream myfile;
    myfile.open (argv[argc-1]);
    myfile << setw(20) << left << "Input_File_Name"
           << setw(15) << "Num_of_Word"
           << setw(20) << "Num_of_Unique_Word"
           << setw(25) << "Num_of_Palindrome_Word"
           << setw(20) << "Num_of_Character"
           << setw(15) << "Num_of_Line" << endl;
    
    for (int i=0; i<argc-2; i++)
    {
        myfile << setw(20) << left << argv[i+1]
               << setw(15) << word_count[i]
               << setw(20) << uniq_word_count[i]
               << setw(25) << palin_word_count[i]
               << setw(20) << char_count[i]
               << setw(15) << line_count[i] << endl;
    }
    myfile.close();
    
    cout << "Job completed! Please find results in " << argv[argc-1]<< endl;

}
