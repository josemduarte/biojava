package org.biojava.nbio.structure.jaligner;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.contact.Pair;
import org.glassfish.jersey.media.multipart.FormDataMultiPart;
import org.glassfish.jersey.media.multipart.MultiPartFeature;
import org.junit.Test;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.Entity;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.io.IOException;

public class TestAlignmentFromTransform {

    private static final String BASE_URL = "https://alignment-west.rcsb.org";
    private static final String SUBMISSION_ENDPOINT = "/api/v1-beta/structures/submit";
    private static final String RESULTS_ENDPOINT = "/api/v1-beta/structures/results?uuid=";

    //-F query='{"mode":"pairwise","method":{"name":"fatcat-rigid"},"structures":[{"entry_id":"1FSL","asym_id":"A"},{"entry_id":"4HHB","asym_id":"A"}]}'

    @Test
    public void testAlignerAgainstRcsbAlignmentService() throws Exception {

        getAlignmentFromRcsbAlignmentService("1FSL", "A", "4HHB", "A");
//        AtomGroupSequence atoms1 = new AtomGroupSequence();
//        AtomGroupSequence atoms2 = new AtomGroupSequence();
//        atoms1.setAtomGroups();
//        atoms2.setAtomGroups();
//        Alignment ali = NeedlemanWunschGotoh.align(atoms1, atoms2, 10, 1);
//        System.out.println(ali.getSequence1());
//        System.out.println(ali.getMarkupLine());
//        System.out.println(ali.getSequence2());
    }

    private Pair<Chain> getAlignmentFromRcsbAlignmentService(String entryId1, String asymId1, String entryId2, String asymId2) throws IOException, InterruptedException {
        String query = String.format(
                "{" +
                    "\"context\":{" +
                    "\"mode\":\"pairwise\"," +
                    "\"method\":{\"name\":\"fatcat-rigid\"}," +
                    "\"structures\":[{\"entry_id\":\"%s\",\"asym_id\":\"%s\"},{\"entry_id\":\"%s\",\"asym_id\":\"%s\"}]" +
                "}}", entryId1, asymId1, entryId2, asymId2);

        FormDataMultiPart multiPart = new FormDataMultiPart();
        multiPart.setMediaType(MediaType.MULTIPART_FORM_DATA_TYPE);

        multiPart.field("query", query);

        Client client = ClientBuilder.newBuilder()
                .register(MultiPartFeature.class)
                .build();
        Response response = client.target(BASE_URL + SUBMISSION_ENDPOINT)
                .request(MediaType.TEXT_PLAIN)
                .header("Content-Type", "multipart/form-data")
                .post(Entity.entity(multiPart, multiPart.getMediaType()));

        if (response.getStatus() !=200)
            throw new IOException("Bad response status: " + response.getStatus());

        String token = response.readEntity(String.class);

        JsonNode node = null;
        long start = System.currentTimeMillis();
        while (System.currentTimeMillis() - start < 30000) {
            response = client.target(BASE_URL + RESULTS_ENDPOINT + token)
                    .request(MediaType.APPLICATION_JSON).get();
            String jsonResp = response.readEntity(String.class);
            ObjectMapper mapper = new ObjectMapper();
            node = mapper.readValue(jsonResp, JsonNode.class);

            String status = node.get("info").get("status").asText();
            if (status.equals("COMPLETE")) {
                break;
            } else {
                Thread.sleep(500);
            }
            node = null;
        }
        if (node == null) throw new IOException("Timed out");
        //System.out.println(node.get("results").toString());
        for (JsonNode oneResult : node.get("results")) {
            System.out.println(oneResult.get("blocks").toString());
        }
        return null;
    }
}
