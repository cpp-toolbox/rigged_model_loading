#ifndef RIGGED_MODEL_LOADING_HPP
#define RIGGED_MODEL_LOADING_HPP

#include <glm/glm.hpp>
#include <string>
#include <utility>
#include <vector>
#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <functional>

#include "sbpt_generated_includes.hpp"

struct VertexBoneData {
    uint indices_of_bones_that_affect_this_vertex[4] = {0};
    float weight_value_of_this_vertex_wrt_bone[4] = {0.0f};

    VertexBoneData() {}

    void add_bone_data(uint BoneID, float Weight);
};

class IVPTRigged {
  public:
    IVPTRigged(std::vector<unsigned int> indices, std::vector<glm::vec3> xyz_positions,
               std::vector<glm::vec2> texture_coordinates, const std::string &texture,
               std::vector<VertexBoneData> bone_data)
        : indices(indices), xyz_positions(xyz_positions), texture_coordinates(texture_coordinates), texture(texture),
          bone_data(bone_data) {};
    Transform transform;
    std::vector<unsigned int> indices;
    std::vector<glm::vec3> xyz_positions;
    std::vector<glm::vec2> texture_coordinates;
    std::string texture;
    std::vector<VertexBoneData> bone_data;
};

// packed version of the above
class IVPTPRigged {
  public:
    IVPTPRigged(std::vector<unsigned int> indices, std::vector<glm::vec3> xyz_positions,
                std::vector<glm::vec2> packed_texture_coordinates, int packed_texture_index, const std::string &texture,
                std::vector<VertexBoneData> bone_data)
        : indices(indices), xyz_positions(xyz_positions), packed_texture_coordinates(packed_texture_coordinates),
          packed_texture_index(packed_texture_index), texture(texture), bone_data(bone_data) {};
    Transform transform;
    std::vector<unsigned int> indices;
    std::vector<glm::vec3> xyz_positions;
    std::vector<glm::vec2> packed_texture_coordinates;
    int packed_texture_index;
    std::string texture;
    std::vector<VertexBoneData> bone_data;
};

struct BoneInfo {
    glm::mat4 local_space_to_bone_space_in_bind_pose_transformation;
    glm::mat4 full_bone_space_to_local_space_transformation = glm::mat4(0);

    BoneInfo(const glm::mat4 &lstbst) { local_space_to_bone_space_in_bind_pose_transformation = lstbst; }
};

class RecIvptRiggedCollector {
  public:
    // we make this here so that we always have access to the importer throughout the lifetime of this object
    // without this there will be segfaults with assimp
    Assimp::Importer importer;
    const aiScene *scene; // for ease
    std::string directory_to_asset_being_loaded;
    std::vector<IVPTRigged> ivptrs;
    void rec_process_nodes(aiNode *node, const aiScene *scene);

    // we work with the transforms of the bones of a mesh, but we never look at the root bone/node's transform
    // if it does have a transform, then in that case it means that the model is already placed somewhere in the world,
    // therefore if we only used the transforms of the non-root nodes it would never get placed where it has to be,
    // therefore since this transform shows how to go from the root node to the origin, then the inverse takes us from
    // the origin to the root node.
    glm::mat4 inverse_root_node_transform;
    int get_next_bone_id(const aiBone *pBone);
    std::unordered_map<std::string, int> bone_name_to_unique_index; // one for each bone.
    std::vector<BoneInfo> bone_unique_idx_to_info;
    std::vector<IVPTRigged> parse_model_into_ivptrs(const std::string &model_path);
    std::vector<VertexBoneData> process_mesh_vertices_bone_data(aiMesh *mesh);
    IVPTRigged process_mesh_ivptrs(aiMesh *mesh, const aiScene *scene);
    // this is used for for eventually binding into uniforms with all the matrices, then in the shader
    // we also have a vertex attribute for each vertex which specifies the id of which matrices to use...
    void set_bone_transforms(float time_in_seconds, std::vector<glm::mat4> &transforms_to_be_set);
    void update_animation_matrices(float animation_time_ticks);
    void rec_update_animation_matrices(float animation_time_ticks, glm::mat4 parent_transform, aiNode *node,
                                       const aiScene *scene);
};

const aiNodeAnim *find_node_anim(const aiAnimation *pAnimation, const std::string &NodeName);

uint find_idx_of_scaling_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
uint find_idx_of_rotation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
uint find_idx_of_translation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_scaling(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_rotation(aiQuaternion &out, float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_translation(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim);

#endif // RIGGED_MODEL_LOADING_HPP
