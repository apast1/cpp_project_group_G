#include<iostream>
#include <string>
using namespace std;

class Person

{
private:
    int age;
    string name;
public:
    void setAge(int);
    void setName(string);
};

class Error
{
public:
    virtual void show()=0;
};
class nameError:public Error
{
public:
    void show()
    {
    cout<<"name is error"<<endl;
    }
 
};

class ageError:public Error
{
public:
    void show()
    {
    cout<<"age is error"<<endl;
    }
};

void Person::setAge(int a)
{
    ageError ag;
    if(a<0||a>100)
        throw ag;
    this->age=a;
}

void Person::setName(string str)
{
    nameError ne;
    if(str=="exit")
        throw ne;
    this->name=str;
}

int main()
{
    Person p;
    
    try
    {
        p.setName("exit");
        p.setAge(101);
    }
    catch(Error &er)
    {
        er.show();
    }
    cout<<"hello world"<<endl;
    
    return 0;
} 
