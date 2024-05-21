#pragma once
#ifndef JOBPARAM_H
#define JOBPARAM_H

/***********************************************************

JobParam 
version 3.4

release note
version3.4:
    ReplaceKey() member is added.
version3.3:
	GetLastError() member is added.
version3.2: (imported from version 2 custom)
	charactors for comment out can be set by user.
version3.1:
    the use of functions _stricmp and _atoi64 are removed.
version3:
	Chainging to the single header file.
version2:
	Functions to set new parameters are added.


************************************************************/


#include <string.h>
#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <locale>



class JobParam{
private:
	std::map<std::string, const char*> m_key_values;
	std::map<std::string, std::vector<const char*> > m_key_multilines;
	
	char* block_begin_head{ nullptr };
	char* block_begin_foot{ nullptr };
	char* block_end_head{ nullptr };
	char* block_end_foot{ nullptr };

	const char* COMMENTOUT = "#\r\n";
	const char* DELIMITER = " \t";



public:
	JobParam(const char* input_file_path, const char* block_begen_word, const char* block_end_word);
	JobParam(const char* input_file_path, const char* block_begen_word, const char* block_end_word, const char* comment_chars);
	JobParam(const char* input_file_path, const char* block_begen_word, const char* block_end_word, const char* comment_chars, const char* comment_word);
	~JobParam();

	int GetInt(const char* key, int default_value) const;
	int GetInt(const char* key) const { return GetInt(key, 0); }
	int64_t GetInt64(const char* key, int64_t default_value) const;
	int64_t GetInt64(const char* key) const {	return GetInt64(key, 0);	}
	double GetDouble(const char* key, double default_value) const;
	double GetDouble(const char* key) const {	return GetDouble(key, 0.0);	}
	int GetIntArray(const char* key, int* buffer, int size) const;
	int GetDoubleArray(const char* key, double* buffer, int size) const;
	const char* GetString(const char* key) const;
	char* DumpString(const char* key) const;
    int IsEqualString(const char* key, const char* compared_str) const;
	//int GetStringArray(const char* key, const char** buffer, int size);
	const char* GetBlockString(const char* key, const int line) const;
	int GetBlockNumLines(const char* key) const;

	int Find(const char* key) const;

	int PrintAll(FILE* fp) const;

    //ver.2機能/////////////////////////////////////////
    int SetInt(const char* key, int value);
    int SetDouble(const char* key, double value);
    int SetIntArray(const char* key, const int* buffer, int size);
    int SetDoubleArray(const char* key, const double* buffer, int size);
    int SetString(const char* key, const char* text);
    int AddBlockString(const char* key, const char* text);
    int Erase(const char* key);

    //ver.3
    bool ReplaceKey(const char* old_key, const char* new_key);

	enum class ErrorCode {
		NO_ERROR, CANNOT_OPEN_FILE, INVALID_BEGIN_END_WORD
	};
private:
	void Open(const char* input_file_path, const char* block_begen_word, const char* block_end_word, const char* comment_chars, const char* comment_word);
	ErrorCode m_last_error = ErrorCode::NO_ERROR;
public:
	ErrorCode GetLastError() {
		return m_last_error;
	};
};


using IteratorKeyValue = std::map<std::string, const char*>::iterator;
using IteratorMultiLine = std::map<std::string, std::vector<const char*> >::iterator;
using IteratorV = std::vector<const char*>::iterator;


//delimiter以外の要素を探す//
template <typename any_char>
any_char* mystrother(any_char* str, const char* delimiter) {
	for (any_char* p = str; *p != '\0'; ++p) {
		const char* r = strchr(delimiter, *p);
		if (!r) {
			return p;
		}
	}
	return NULL;
}

//delimiterに含まれる要素を探す//
template <typename any_char>
any_char* mystrany(any_char* str, const char* delimiter) {
	for (any_char* p = str; *p != '\0'; ++p) {
		const char* r = strchr(delimiter, *p);
		if (r) {
			return p;
		}
	}
	return NULL;
}

template <typename any_char>
char* mystrdup(any_char* str) {
	size_t length = strlen(str);
	if (length > 0) {
		char* buffer = new char[length + 1];
		strcpy(buffer, str);
		return buffer;
	} else {
		return NULL;
	}
}


//末尾からみてdelimiterに含まれる最後の要素をトリム//
template <typename any_char>
any_char* mystrrtrim(any_char* str, const char* delimiter) {

	for (any_char* p = str + strlen(str) - 1; p != str - 1; p--) {
		const char* r = strchr(delimiter, *p);
		if (r) {
			*p = '\0';
		} else {
			return p;
		}
	}
	return NULL;
}


inline int mystrrcmp(const char* target, const char* word) {
	if (strlen(target) >= strlen(word)) {		
		return strcmp(target + strlen(target) - strlen(word), word);
	}

	return -1;

}

inline JobParam::JobParam(const char* input_file_path, const char* block_begen_word, const char* block_end_word) {
	Open(input_file_path, block_begen_word, block_end_word, COMMENTOUT, nullptr);
}

inline JobParam::JobParam(const char* input_file_path, const char* block_begen_word, const char* block_end_word, const char* comment_chars){
	Open(input_file_path, block_begen_word, block_end_word, comment_chars, nullptr);
}

inline JobParam::JobParam(const char* input_file_path, const char* block_begen_word, const char* block_end_word, const char* comment_chars, const char* comment_word) {
	Open(input_file_path, block_begen_word, block_end_word, comment_chars, comment_word);
}

inline
void JobParam::Open(const char* input_file_path, const char* block_begen_word, const char* block_end_word, const char* comment_chars, const char* comment_word) {


	FILE* fp = fopen(input_file_path, "r");
	if (!fp) {
		m_last_error = ErrorCode::CANNOT_OPEN_FILE;
		return;
	}

	//ブロックステートメントの判定ワード抽出//
	{
		const char* p = strchr(block_begen_word, '*');
		if (p) {
			if (p != block_begen_word) {
				size_t length = p - block_begen_word;
				block_begin_head = new char[length + 1];
				strncpy(block_begin_head, block_begen_word, length);
				block_begin_head[length] = '\0';
			} else {
				block_begin_head = new char[1];
				*block_begin_head = '\0';
			}
			p++;
			if (*p != '\0') {
				size_t length = strlen(p);
				block_begin_foot = new char[length + 1];
				strcpy(block_begin_foot, p);
			} else {
				block_begin_foot = new char[1];
				*block_begin_foot = '\0';
			}
		} else {
			fclose(fp);
			m_last_error = ErrorCode::INVALID_BEGIN_END_WORD;
			return;
		}
	}

	{
		const char* p = strchr(block_end_word, '*');
		if (p) {
			if (p != block_end_word) {
				size_t length = p - block_end_word;
				block_end_head = new char[length + 1];
				strncpy(block_end_head, block_end_word, length);
				block_end_head[length] = '\0';
			} else {
				block_end_head = new char[1];
				*block_end_head = '\0';
			}
			p++;
			if (*p != '\0') {
				size_t length = strlen(p);
				block_end_foot = new char[length + 1];
				strcpy(block_end_foot, p);
			} else {
				block_end_foot = new char[1];
				*block_end_foot = '\0';
			}
		} else {
			fclose(fp);
			m_last_error = ErrorCode::INVALID_BEGIN_END_WORD;
			return;
		}
	}


	const int SIZE_ROW = 1024;
	char line[SIZE_ROW];

	bool is_in_comment = false;

	while (fgets(line, SIZE_ROW, fp) != NULL) {

		//複数行コメントアウト//
		if (is_in_comment) {
			if (strncmp(line, "*/", 2) == 0) {
				is_in_comment = false;
			}
			continue;
		} else {
			if (strncmp(line, "/*", 2) == 0) {
				is_in_comment = true;
				continue;
			}
		}

		if(comment_word){//コメントを削除, 単語単位("//"や"<VPS>"など), 右側トリム//
			char* p = strstr(line, comment_word);
			if (p) { *p = '\0'; }
			mystrrtrim(line, DELIMITER);
		}

		{//コメントを削除, 文字単位, 右側トリム//
			char* p = mystrany(line, comment_chars);
			if (p) { *p = '\0'; }
			mystrrtrim(line, DELIMITER);
		}

		char* key_head = mystrother(line, DELIMITER);
		if (key_head) {

			char* value_head = NULL;

			char* key_end = mystrany(key_head + 1, DELIMITER);
			if (key_end) {
				/*int length = (int)(key_end - key_head);
				key = new char[length + 1];
				strncpy(key, key_head, length);
				key[length] = '\0';
				*/
				*key_end = '\0';
				value_head = mystrother(key_end + 1, DELIMITER);


			} else {
				/*
				int length = strlen(key_head);
				key = new char[length + 1];
				strcpy(key, key_head);
				*/
				//keyだけ見つかってバリューが無い//
			}

			//Begin-Endで挟まるマルチ行指定かどうか
			int flag_multiline = 1;
			if (strlen(block_begin_head)>0) {
				if (strncmp(key_head, block_begin_head, strlen(block_begin_head)) != 0) {
					flag_multiline = 0;
				}
			}
			if (strlen(block_begin_foot)>0) {
				if (mystrrcmp(key_head, block_begin_foot) != 0) {
					flag_multiline = 0;
				}
			}

			//マルチ行でない//
			if (flag_multiline == 0) {
				//No: 通常のkeyとvalueを追加//
				if (value_head) {
					m_key_values.insert(std::pair<std::string, const char*>(key_head, mystrdup(value_head)));
				}

			} else {//Yes: .Endまでの複数行を追加//

					//※これ以降はlineが書き変えられるためkey_headの値を保持されない//


					//末行の文字列を作成
				size_t length = strlen(key_head) - strlen(block_begin_head) - strlen(block_begin_foot);
				char* key = new char[length + 1];
				char* term_key = new char[length + strlen(block_end_head) + strlen(block_end_foot) + 1];
				strncpy(key, key_head + strlen(block_begin_head), length);
				key[length] = '\0';
				if (strlen(block_end_head) > 0) { strcpy(term_key, block_end_head); }
				strcpy(term_key + strlen(block_end_head), key);
				if (strlen(block_end_foot) > 0) { strcpy(term_key + strlen(term_key), block_end_foot); }


				//追加するvectorクラス//
				std::vector<const char*> multilines;


				//新たに行を読み込む

				while (fgets(line, SIZE_ROW, fp) != NULL) {

					//複数行コメントアウト//
					if (is_in_comment) {
						if (strncmp(line, "*/", 2) == 0) {
							is_in_comment = false;
						}
						continue;
					} else {
						if (strncmp(line, "/*", 2) == 0) {
							is_in_comment = true;
							continue;
						}
					}

					{//コメントを削除, 右側トリム//
						char* p = mystrany(line, comment_chars);
						if (p) { *p = '\0'; }
						mystrrtrim(line, DELIMITER);
					}

					char* line_head = mystrother(line, DELIMITER);
					if (line_head) {
						if (strcmp(line_head, term_key) == 0) {
							//末行を見つけた//

							m_key_multilines.insert(std::pair<std::string, std::vector<const char*> >(key, multilines));

							break;

						} else {
							//末行ではなかったのでバッファに追加//
							multilines.push_back(mystrdup(line_head));
						}
					}

				}


				delete[] term_key;
				delete[] key;


			}//Yes: .Endまでの複数行を追加//
		}//key_head//


	}

	fclose(fp);

	
}


inline JobParam::~JobParam() {

	//全ての要素を削除//
	for (IteratorKeyValue it = m_key_values.begin(); it != m_key_values.end(); it++) {
		delete[](it->second);
	}
	m_key_values.clear();

	//全てのマルチ行要素を削除//
	for (IteratorMultiLine it = m_key_multilines.begin(); it != m_key_multilines.end(); it++) {
		for (IteratorV k = it->second.begin(); k != it->second.end(); k++) {
			delete[] * k;
		}
		it->second.clear();
	}
	m_key_multilines.clear();


	if (block_begin_head) { delete[] block_begin_head; }
	if (block_begin_foot) { delete[] block_begin_foot; }
	if (block_end_head) { delete[] block_end_head; }
	if (block_end_foot) { delete[] block_end_foot; }

}


inline int JobParam::GetInt(const char* key, int default_value) const {

	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return default_value;
	} else {
		// キーが "Key" の要素がある場合の処理
		return atoi(it->second);
	}

}

inline int64_t JobParam::GetInt64(const char* key, int64_t default_value) const {

	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return default_value;
	} else {
		// キーが "Key" の要素がある場合の処理
		return (int64_t)strtol(it->second,nullptr,10);
	}

}

inline int JobParam::GetIntArray(const char* key, int* buffer, int size) const {
	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return 0;
	} else {
		// キーが "Key" の要素がある場合の処理
		char* temp = mystrdup(it->second);
		int count = 0;
		char* tp = strtok(temp, DELIMITER);
		while (tp) {
			buffer[count] = atoi(tp);
			count++;
			if (count >= size) break;
			tp = strtok(NULL, DELIMITER);
		}
		delete[] temp;
		return count;
	}

}

inline double JobParam::GetDouble(const char* key, double default_value) const {

	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return default_value;
	} else {
		// キーが "Key" の要素がある場合の処理
		return strtod(it->second, NULL);
	}

}

inline int JobParam::GetDoubleArray(const char* key, double* buffer, int size) const {
	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return 0;
	} else {
		// キーが "Key" の要素がある場合の処理
		char* temp = mystrdup(it->second);
		int count = 0;
		char* tp = strtok(temp, DELIMITER);
		while (tp) {
			buffer[count] = strtod(tp, NULL);
			count++;
			if (count >= size) break;
			tp = strtok(NULL, DELIMITER);
		}
		delete[] temp;
		return count;
	}

}


inline const char* JobParam::GetString(const char* key) const {

	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return NULL;
	} else {
		// キーが "Key" の要素がある場合の処理
		return it->second;
	}

}

inline char* JobParam::DumpString(const char* key) const {
	auto it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		return NULL;
	} else {
		// キーが "Key" の要素がある場合の処理
		return mystrdup(it->second);
	}

}

/*
引数の文字列と、valueの文字列が同じかどうか比較
同じなら1を返す。
大文字小文字は区別しない
*/
inline int JobParam::IsEqualString(const char* key, const char* compared_str) const {
	const char* value = GetString(key);
	if (value && compared_str) {
		std::string vlower(value);
		for (auto&& c : vlower) {
			c = std::tolower(c);
		}

		std::string cslower(compared_str);
		for (auto&& c : cslower) {
			c = std::tolower(c);
		}

		if (vlower == cslower) {
			return 1;
		}
	}

	return 0;

}

/*
beginとendで囲まれたブロック要素の中身を、一行ずつ取得する。
line = 0,1,2...で指定した行の文字列を返す
lineで指定した行がブロック要素を超えていたらNULLを返す
*/
inline const char* JobParam::GetBlockString(const char* key, const int line) const {

	auto it = m_key_multilines.find(key);
	if (it == m_key_multilines.end()) {
		return NULL;
	} else {
		// キーが "Key" の要素がある場合の処理
		if (line >= it->second.size()) {
			return NULL;
		} else {
			return it->second[line];
		}


	}

}

/*
beginとendで囲まれたブロック要素の行数を返す
*/
inline int JobParam::GetBlockNumLines(const char* key) const {

	auto it = m_key_multilines.find(key);
	if (it == m_key_multilines.end()) {
		return 0;
	} else {
		return (int)it->second.size();

	}
}

inline int JobParam::Find(const char* key) const {

	auto it = m_key_values.find(key);
	if (m_key_values.find(key) != m_key_values.end()) {
		return 1;
	} else if (m_key_multilines.find(key) != m_key_multilines.end()) {
		return 2;
	} else {
		return 0;
	}
}

inline int JobParam::PrintAll(FILE* fp) const {
	if (!fp) return -1;

	int count = 0;

	//全ての要素を削除//
	for (auto it = m_key_values.begin(); it != m_key_values.end(); it++) {
		fprintf(fp, "%s\t%s\n", it->first.c_str(), it->second);
		count++;
	}
	//m_key_values.clear();

	//全てのマルチ行要素を削除//
	for (auto it = m_key_multilines.begin(); it != m_key_multilines.end(); it++) {
		fprintf(fp, "%s%s%s\n", block_begin_head, it->first.c_str(), block_begin_foot);
		for (auto k = it->second.begin(); k != it->second.end(); k++) {
			fprintf(fp, "\t%s\n", *k);
			count++;
		}
		fprintf(fp, "%s%s%s\n", block_end_head, it->first.c_str(), block_end_foot);
		count += 2;
	}

	return count;

}


//ver.2機能/////////////////////////////////////////////////

inline int JobParam::SetInt(const char* key, int value) {

	char text[16];
	sprintf(text, "%d", value);

	return SetString(key, text);

}

inline int JobParam::SetIntArray(const char* key, const int* buffer, int size) {

	char* text = new char[16 * size];
	for (int i = 0; i < size; i++) {
		sprintf(text + strlen(text), " %d", buffer[i]);
	}

	return SetString(key, text);

}

inline int JobParam::SetDouble(const char* key, double value) {

	char text[16];
	sprintf(text, "%lf", value);

	return SetString(key, text);

}

inline int JobParam::SetDoubleArray(const char* key, const double* buffer, int size) {
	char* text = new char[16 * size];
	for (int i = 0; i < size; i++) {
		sprintf(text + strlen(text), " %lf", buffer[i]);
	}

	return SetString(key, text);

}


inline int JobParam::SetString(const char* key, const char* text) {

	IteratorKeyValue it = m_key_values.find(key);
	if (it == m_key_values.end()) {
		m_key_values.insert(std::pair<std::string, const char*>(key, mystrdup(text)));
		return 1;
	} else {
		// キーが "Key" の要素がある場合の処理
		delete[] it->second;
		it->second = mystrdup(text);
		return 1;
	}

}


inline int JobParam::AddBlockString(const char* key, const char* text) {

	IteratorMultiLine it = m_key_multilines.find(key);
	if (it == m_key_multilines.end()) {

		std::vector<const char*> multilines;

		multilines.push_back(mystrdup(text));

		m_key_multilines.insert(std::pair<std::string, std::vector<const char*> >(key, multilines));
		return 1;
	} else {
		// キーが "Key" の要素がある場合の処理
		it->second.push_back(mystrdup(text));
		return 1;
	}

}

inline int JobParam::Erase(const char* key) {

	{
		IteratorKeyValue it = m_key_values.find(key);
		if (it != m_key_values.end()) {
			delete[] it->second;
			m_key_values.erase(it);
			return 1;
		};
	}

	{
		IteratorMultiLine it = m_key_multilines.find(key);
		if (it != m_key_multilines.end()) {
			for (int i = 0; i < it->second.size(); i++) {
				delete[] it->second[i];
			}
			m_key_multilines.erase(it);
			return 1;
		};
	}

	return 0;


}


inline bool JobParam::ReplaceKey(const char* old_key, const char* new_key) {

    {
        IteratorKeyValue it = m_key_values.find(old_key);
        if (it != m_key_values.end()) {
            m_key_values.insert(std::make_pair(new_key, it->second));
            m_key_values.erase(it);
            return true;
        };
    }

    {
        IteratorMultiLine it = m_key_multilines.find(old_key);
        if (it != m_key_multilines.end()) {
            m_key_multilines.insert(std::make_pair(new_key, it->second));
            m_key_multilines.erase(it);
            return true;
        };
    }
    return false;
}



#endif	//!JOBPARAM_H
