#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

struct vec3 {
public:
    double coord[3];
    double x() const {return this->coord[0];}
    double y() const {return this->coord[1];}
    double z() const {return this->coord[2];}
    vec3(double x, double y, double z) {
        this->coord[0] = x;
        this->coord[1] = y;
        this->coord[2] = z;
    }
    vec3 cross(const vec3& vec) const {
        vec3 mult = vec3(this->coord[2] * vec.coord[1] - this->coord[1] * vec.coord[2],
                         this->coord[0] * vec.coord[2] - this->coord[2] * vec.coord[0],
                         this->coord[1] * vec.coord[0] - this->coord[0] * vec.coord[1]);
        return mult;
    }
    vec3& operator+=(const vec3& other) {
        this->coord[0] += other.coord[0];
        this->coord[1] += other.coord[1];
        this->coord[2] += other.coord[2];
        return *this;
    }
    vec3 operator+(const vec3& vec) const {
        return vec3(this->x() + vec.x(), this->y() + vec.y(), this->z() + vec.z());
    }
    vec3 operator-(const vec3& vec) const {
        return vec3(this->x() - vec.x(), this->y() - vec.y(), this->z() - vec.z());
    }
    double dot(const vec3& vec) const {
        return this->coord[0] * vec.coord[0] + this->coord[1] * vec.coord[1] + this->coord[2] * vec.coord[2];
    }
    double length() const {
        double sq_len = this->dot(*this);
        return std::sqrt(sq_len);
    }
    vec3 vector_projection(const vec3&) const;
};

vec3 operator*(double t, const vec3& vec) {
    return vec3(t * vec.x(), t * vec.y(), t * vec.z());
}

vec3 operator/(const vec3& vec, double t) {
    return vec3(vec.x() / t, vec.y() / t, vec.z() / t);
}

vec3 unit_vector(const vec3& vec) {
    return vec / vec.length();
}

vec3 vec3::vector_projection(const vec3& vec) const {
    if (vec.length() < 10e-6) return vec3(0, 0, 0);
    return this->dot(vec) * unit_vector(*this);
}

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    vec3 pos;
    // Некая величина, в попугаях
    double smth;
    // Скорость
    vec3 vel;

public:
    // Конструктор по умолчанию
    CalcNode() : pos(vec3(0, 0, 0)), smth(smth), vel(vec3(0, 0, 0))
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
            : pos(vec3(x, y, z)), smth(smth), vel(vec3(vx, vy, vz))
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau, vec3& w_cap, vec3& w_low) {
        vel = w_cap.cross(pos - vec3(0, 5, 1)) + w_low.cross(pos - vec3(0, 5, 1));
        pos += tau * vel;
    }

    void heat(double time, vec3& w_cap, vec3& w_low) {
        smth = std::sin((pos - pos.vector_projection(w_cap) - pos.vector_projection(w_low)).length() - 30 * time);
    }
};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            // Модельная скалярная величина распределена как-то вот так
            double smth = std::sin(pointY + pointZ);
            nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(int step, double tau, vec3& w_cap, vec3& w_low) {
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].move(tau, w_cap, w_low);
            nodes[i].heat(step * tau, w_cap, w_low);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].pos.coord[0], nodes[i].pos.coord[1], nodes[i].pos.coord[2]);

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {nodes[i].vel.coord[0], nodes[i].vel.coord[1], nodes[i].vel.coord[2]};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            smth->InsertNextValue(nodes[i].smth);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "./res/seal-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    // Шаг точек по пространству
    double h = 4.0;
    // Шаг по времени
    double tau = 0.001;

    const unsigned int GMSH_TETR_CODE = 4;

    // Теперь придётся немного упороться:
    // (а) построением сетки средствами gmsh,
    // (б) извлечением данных этой сетки в свой код.
    gmsh::initialize();
    gmsh::model::add("seal");

    // Считаем STL
    try {
        gmsh::merge("./stl/Popo.stl"); 
        // путь к файлу отсчитывается от точки запуска 
        // если вы собирали все в директории build как цивилизованные люди, 
        // то перейдите на уровень выше
        // и запускайте бинарник как ./build/tetr3d
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Восстановим геометрию
    double angle = 120;
    bool forceParametrizablePatches = true;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.3);

    // Зададим объём по считанной поверхности
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    // Зададим мелкость желаемой сетки
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    // Построим сетку
    gmsh::model::mesh::generate(3);

    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();
    
    double freq_mag_cap = 3.0;
    double freq_round_low = 3.0;
    double freq_mag_low = 20.0;
    double max_mag_cap = 1.0, cur_mag_low = 0.0;
    vec3 w_cap = vec3(0, 0, max_mag_cap);
    vec3 w_low = vec3(0, 0, 0);
    vec3 unit_w_low = vec3(1, 0, 0);
    for(int step = 1; step <= 5000; step++) {
        w_low = cur_mag_low * unit_w_low;
        mesh.doTimeStep(step, tau, w_cap, w_low);
        if (step % 10 == 0) mesh.snapshot(step);
        w_cap.coord[2] -= max_mag_cap * freq_mag_cap * std::sin(freq_mag_cap * step * tau) * tau;
        cur_mag_low += freq_round_low * freq_mag_low * std::cos(freq_mag_low * step * tau) * tau;
        unit_w_low += tau * w_cap.cross(unit_w_low);
        if (step % 100 == 0) cout << step << endl;
    }

    return 0;
}
